#!/usr/bin/env Rscript
# pseudobulk.R — DESeq2 pseudo-bulk differential expression wrapper for Nextflow DSL2
#
# Usage:
#   Rscript pseudobulk.R \
#     --input      <clustered.h5ad> \
#     --groupby    <leiden> \
#     --sample-key <sample> \
#     --min-cells  10 \
#     --output     <de_results.csv>
#
# Algorithm: one-vs-rest pseudo-bulk DESeq2 per cluster.
#   Aggregation: sum raw counts (counts assay) per (cluster × sample).
#   Design: ~ group  (group = 1 if sample belongs to cluster, else 0).
#   Minimum: ≥2 samples with ≥min_cells cells in either arm; cluster skipped otherwise.
#   Test: Wald test for group1 vs group0 coefficient.
#   Score: Wald statistic (stat column from DESeq2 results).
#   p-values corrected with BH across all genes for each cluster.
#
# Output columns: group, gene, score, logfoldchange, pval, pval_adj
#
# Requires: bioconductor-deseq2 (conda: bioconductor-deseq2 >= 1.40)

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(DESeq2)
  library(Matrix)
})

parse_args_pseudobulk <- function() {
  option_list <- list(
    make_option("--input",      type = "character",
                help = "Input .h5ad path (must have assay 'counts')"),
    make_option("--groupby",    type = "character", default = "leiden",
                help = "obs column with cluster labels [default: leiden]"),
    make_option("--sample-key", type = "character", default = "sample",
                help = "obs column with sample IDs [default: sample]"),
    make_option("--min-cells",  type = "integer",   default = 10L,
                help = "Min cells per (cluster, sample) pseudo-bulk [default: 10]"),
    make_option("--output",     type = "character",
                help = "Output CSV path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$input) || is.null(opts$output)) {
    print_help(parser)
    stop("--input and --output are required", call. = FALSE)
  }
  opts
}

# Aggregate raw counts into pseudo-bulk matrix: genes x (cluster__sample) pseudoreplicates.
# Returns list: mat (genes x samples matrix), meta (data.frame with cluster and sample cols).
aggregate_pseudobulk <- function(counts_mat, clusters, samples) {
  groups <- paste(clusters, samples, sep = "__")
  unique_groups <- sort(unique(groups))

  agg <- vapply(unique_groups, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1L) {
      as.numeric(counts_mat[, idx, drop = FALSE])
    } else {
      Matrix::rowSums(counts_mat[, idx, drop = FALSE])
    }
  }, numeric(nrow(counts_mat)))

  rownames(agg) <- rownames(counts_mat)
  colnames(agg) <- unique_groups

  parts <- strsplit(unique_groups, "__", fixed = TRUE)
  meta  <- data.frame(
    cluster    = vapply(parts, `[[`, character(1L), 1L),
    sample     = vapply(parts, function(x) paste(x[-1], collapse = "__"), character(1L)),
    row.names  = unique_groups,
    stringsAsFactors = FALSE
  )
  list(mat = agg, meta = meta)
}

run_pseudobulk <- function(input_path, groupby, sample_key, min_cells, output_path) {
  message("pseudobulk.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  # Validate inputs
  for (col in c(groupby, sample_key)) {
    if (!col %in% names(colData(sce))) {
      stop("colData missing '", col, "'. Available: ",
           paste(names(colData(sce)), collapse = ", "))
    }
  }
  if (!"counts" %in% assayNames(sce)) {
    stop("Input .h5ad must contain assay 'counts'")
  }

  counts_mat  <- assay(sce, "counts")           # genes x cells (sparse ok)
  clusters    <- as.character(colData(sce)[[groupby]])
  samples     <- as.character(colData(sce)[[sample_key]])
  cluster_ids <- sort(unique(clusters))

  # Build full pseudo-bulk aggregation once
  pb        <- aggregate_pseudobulk(counts_mat, clusters, samples)
  pb_mat    <- pb$mat
  pb_meta   <- pb$meta

  # Cell counts per (cluster, sample) for filtering
  cell_counts <- table(clusters, samples)

  all_results <- vector("list", length(cluster_ids))

  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    message("pseudobulk.R: testing cluster ", cl,
            " (", i, "/", length(cluster_ids), ")")

    # Assign group label to each pseudoreplicate
    pb_meta$group <- ifelse(pb_meta$cluster == cl, "1", "0")

    # Identify valid samples: ≥min_cells cells in this cluster or rest
    in_samples  <- rownames(cell_counts)[cell_counts[cl, ] >= min_cells]
    rest_mask   <- pb_meta$group == "0"
    out_samples <- pb_meta$sample[rest_mask &
                     vapply(pb_meta$sample[rest_mask], function(s) {
                       total_s <- sum(cell_counts[, s])
                       total_s >= min_cells
                     }, logical(1L))]
    out_samples <- unique(out_samples)

    n_in  <- length(in_samples)
    n_out <- length(out_samples)
    if (n_in < 2L || n_out < 2L) {
      message("pseudobulk.R:  cluster ", cl, " skipped (",
              n_in, " in-group samples, ", n_out, " out-group samples)")
      next
    }

    keep_cols <- rownames(pb_meta)[pb_meta$group == "1" & pb_meta$sample %in% in_samples  |
                                   pb_meta$group == "0" & pb_meta$sample %in% out_samples]

    sub_mat  <- round(pb_mat[, keep_cols, drop = FALSE])
    sub_meta <- pb_meta[keep_cols, , drop = FALSE]
    sub_meta$group <- factor(sub_meta$group, levels = c("0", "1"))

    # Remove all-zero genes for this subset
    keep_genes <- rowSums(sub_mat) > 0
    if (sum(keep_genes) == 0L) {
      message("pseudobulk.R:  cluster ", cl, " skipped (all genes zero)")
      next
    }
    sub_mat <- sub_mat[keep_genes, , drop = FALSE]

    dds <- tryCatch({
      dds_obj <- DESeqDataSetFromMatrix(
        countData = sub_mat,
        colData   = sub_meta,
        design    = ~ group
      )
      suppressMessages(DESeq(dds_obj, quiet = TRUE))
    }, error = function(e) {
      message("pseudobulk.R:  DESeq2 failed for cluster ", cl, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(dds)) next

    res <- tryCatch(
      as.data.frame(results(dds, name = "group_1_vs_0", independentFiltering = FALSE)),
      error = function(e) {
        message("pseudobulk.R:  results() failed for cluster ", cl, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res)) next

    # Fill NA stats with 0 / 1
    res$stat[is.na(res$stat)]   <- 0
    res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
    res$pvalue[is.na(res$pvalue)] <- 1

    all_results[[i]] <- data.frame(
      group         = cl,
      gene          = rownames(res),
      score         = res$stat,
      logfoldchange = res$log2FoldChange,
      pval          = res$pvalue,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, Filter(Negate(is.null), all_results))
  if (is.null(df) || nrow(df) == 0L) {
    stop("pseudobulk.R: no clusters produced results.")
  }

  # BH correction per cluster
  df$pval_adj <- ave(df$pval, df$group,
                     FUN = function(p) p.adjust(p, method = "BH"))

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    df[, c("group", "gene", "score", "logfoldchange", "pval", "pval_adj")],
    output_path, row.names = FALSE, quote = FALSE
  )
  message("pseudobulk.R: wrote ", nrow(df), " rows -> ", output_path)
}

main <- function() {
  opts <- parse_args_pseudobulk()
  run_pseudobulk(opts$input, opts$groupby, opts[["sample-key"]],
                 opts[["min-cells"]], opts$output)
}

if (!interactive()) main()
