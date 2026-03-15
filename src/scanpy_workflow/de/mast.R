#!/usr/bin/env Rscript
# mast.R — MAST hurdle-model differential expression wrapper for Nextflow DSL2
#
# Usage:
#   Rscript mast.R \
#     --input   <clustered.h5ad> \
#     --groupby <leiden> \
#     --output  <de_results.csv>
#
# Algorithm: one-vs-rest hurdle model per cluster.
#   Expression: log_norm assay (log1p-normalised, genes x cells).
#   Covariates: cellular detection rate (ngeneson).
#   Test: likelihood ratio test for the cluster assignment coefficient.
#   Score: sign(LFC) * sqrt(lambda_hurdle)  — signed chi-squared root.
#   p-values corrected with BH across all genes for each cluster.
#
# Output columns: group, gene, score, logfoldchange, pval, pval_adj
#
# Requires: bioconductor-mast (conda: bioconductor-mast >= 1.26)

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(MAST)
  library(Matrix)
})

parse_args_mast <- function() {
  option_list <- list(
    make_option("--input",   type = "character",
                help = "Input .h5ad path (must have assay 'log_norm')"),
    make_option("--groupby", type = "character", default = "leiden",
                help = "obs column with cluster labels [default: leiden]"),
    make_option("--output",  type = "character",
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

run_mast <- function(input_path, groupby, output_path) {
  message("mast.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!groupby %in% names(colData(sce))) {
    stop("colData missing '", groupby, "'. Available: ",
         paste(names(colData(sce)), collapse = ", "))
  }
  if (!"log_norm" %in% assayNames(sce)) {
    stop("Input .h5ad must contain assay 'log_norm'")
  }

  expr_mat    <- as.matrix(assay(sce, "log_norm"))   # genes x cells
  clusters    <- as.character(colData(sce)[[groupby]])
  cluster_ids <- sort(unique(clusters))

  # Cellular detection rate: fraction of genes expressed per cell
  ngeneson <- colMeans(expr_mat > 0)

  all_results <- vector("list", length(cluster_ids))

  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    message("mast.R: testing cluster ", cl,
            " (", i, "/", length(cluster_ids), ")")

    group_vec <- as.integer(clusters == cl)   # 1 = in cluster, 0 = rest

    cdat <- data.frame(
      wellKey  = colnames(sce),
      group    = factor(group_vec, levels = c(0L, 1L)),
      ngeneson = ngeneson,
      row.names = colnames(sce),
      stringsAsFactors = FALSE
    )
    fdat <- data.frame(
      primerid  = rownames(sce),
      row.names = rownames(sce)
    )

    sca <- suppressMessages(
      FromMatrix(expr_mat, cData = cdat, fData = fdat)
    )

    # Fit hurdle model: ~ group + ngeneson
    zlm_fit <- tryCatch(
      suppressMessages(zlm(~ group + ngeneson, sca)),
      error = function(e) {
        message("mast.R:  zlm failed for cluster ", cl, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(zlm_fit)) next

    # LRT for group coefficient
    lrt <- tryCatch(
      suppressMessages(lrTest(zlm_fit, "group")),
      error = function(e) {
        message("mast.R:  lrTest failed for cluster ", cl, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(lrt)) next

    # Hurdle p-value (combined continuous + discrete components)
    hurdle_pval <- lrt[, "hurdle", "Pr(>Chisq)"]
    hurdle_pval[is.na(hurdle_pval)] <- 1

    # Log fold change: mean log_norm in cluster vs rest
    in_cl <- clusters == cl
    lfc   <- rowMeans(expr_mat[,  in_cl, drop = FALSE]) -
             rowMeans(expr_mat[, !in_cl, drop = FALSE])

    # Score: signed square root of LRT statistic (sign given by LFC direction)
    lambda <- lrt[, "hurdle", "lambda"]
    lambda[is.na(lambda)] <- 0
    score  <- sign(lfc) * sqrt(pmax(lambda, 0))

    all_results[[i]] <- data.frame(
      group         = cl,
      gene          = rownames(sce),
      score         = score,
      logfoldchange = lfc,
      pval          = hurdle_pval,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, Filter(Negate(is.null), all_results))
  if (is.null(df) || nrow(df) == 0) {
    stop("mast.R: no clusters produced results.")
  }

  # BH correction per cluster
  df$pval_adj <- ave(df$pval, df$group,
                     FUN = function(p) p.adjust(p, method = "BH"))

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    df[, c("group", "gene", "score", "logfoldchange", "pval", "pval_adj")],
    output_path, row.names = FALSE, quote = FALSE
  )
  message("mast.R: wrote ", nrow(df), " rows -> ", output_path)
}

main <- function() {
  opts <- parse_args_mast()
  run_mast(opts$input, opts$groupby, opts$output)
}

if (!interactive()) main()
