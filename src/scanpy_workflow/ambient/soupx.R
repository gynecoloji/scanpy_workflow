#!/usr/bin/env Rscript
# soupx.R — SoupX ambient RNA removal wrapper for Nextflow DSL2
#
# Usage:
#   Rscript soupx.R \
#     --filtered-h5ad <filtered.h5ad> \
#     --raw-dir       <cellranger_raw_feature_bc_matrix/> \
#     --output        <decontaminated.h5ad>
#
# Reads the raw CellRanger directory with a base R helper (no Seurat required).
# Writes corrected integer counts + recomputed norm/log_norm/X layers.

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SoupX)
  library(Matrix)
})

parse_args_soupx <- function() {
  option_list <- list(
    make_option("--filtered-h5ad", type = "character", dest = "filtered_h5ad",
                help = "Filtered .h5ad input path"),
    make_option("--raw-dir", type = "character", dest = "raw_dir",
                help = "CellRanger raw_feature_bc_matrix directory"),
    make_option("--output", type = "character", help = "Output .h5ad path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$filtered_h5ad) || is.null(opts$raw_dir) || is.null(opts$output)) {
    print_help(parser)
    stop("--filtered-h5ad, --raw-dir, and --output are required", call. = FALSE)
  }
  opts
}

#' Read a 10x CellRanger matrix directory (gzipped files) without Seurat.
#' Returns a genes x barcodes sparse matrix (dgCMatrix).
read_10x_matrix <- function(dir) {
  bc_con   <- gzfile(file.path(dir, "barcodes.tsv.gz"), "r")
  barcodes <- readLines(bc_con)
  close(bc_con)

  ft_con     <- gzfile(file.path(dir, "features.tsv.gz"), "r")
  feat_lines <- readLines(ft_con)
  close(ft_con)
  feature_ids <- sub("\t.*", "", feat_lines)

  mtx_con <- gzfile(file.path(dir, "matrix.mtx.gz"), "r")
  mtx     <- Matrix::readMM(mtx_con)
  close(mtx_con)

  rownames(mtx) <- feature_ids
  colnames(mtx) <- barcodes
  as(mtx, "dgCMatrix")
}

run_soupx <- function(filtered_h5ad, raw_dir, output_path) {
  message("soupx.R: reading filtered .h5ad from ", filtered_h5ad)
  sce_filt <- readH5AD(filtered_h5ad, reader = "R")

  if (!"counts" %in% assayNames(sce_filt)) {
    stop("Filtered .h5ad must contain assay 'counts'")
  }
  toc <- assay(sce_filt, "counts")   # table of counts (filtered cells)

  message("soupx.R: reading raw CellRanger directory from ", raw_dir)
  tod <- read_10x_matrix(raw_dir)    # table of droplets (all barcodes incl. empty)

  # Align features: raw dir may have a different feature set
  common_genes <- intersect(rownames(tod), rownames(toc))
  if (length(common_genes) == 0) {
    stop("No common genes between raw dir and filtered .h5ad. Check gene name format.")
  }
  tod <- tod[common_genes, ]
  toc <- toc[common_genes, ]

  # Verify that raw matrix contains all filtered barcodes
  missing_bc <- setdiff(colnames(toc), colnames(tod))
  if (length(missing_bc) > 0) {
    stop(sprintf(
      "soupx.R: %d filtered barcodes not found in raw dir (e.g. %s). Check that --raw-dir points to the correct CellRanger raw_feature_bc_matrix.",
      length(missing_bc), paste(head(missing_bc, 3), collapse = ", ")
    ))
  }

  message("soupx.R: building SoupChannel and estimating contamination")
  sc <- SoupChannel(tod, toc)

  # Use leiden clusters if available; otherwise quick-cluster from lib sizes
  if ("leiden" %in% names(colData(sce_filt))) {
    cl <- setNames(as.character(colData(sce_filt)[colnames(toc), "leiden"]),
                   colnames(toc))
  } else {
    lib <- colSums(toc)
    cl  <- setNames(ifelse(lib > median(lib), "1", "0"), colnames(toc))
  }
  sc <- setClusters(sc, cl)
  sc <- autoEstCont(sc, doPlot = FALSE)

  message("soupx.R: adjusting counts")
  corrected <- adjustCounts(sc, roundToInt = TRUE)
  # adjustCounts may change column order; re-align to toc column order
  corrected <- corrected[, colnames(toc)]

  # Recompute downstream layers from corrected counts
  lib_sizes    <- colSums(corrected)
  lib_sizes[lib_sizes == 0] <- 1
  norm_mat     <- sweep(corrected, 2, lib_sizes / 1e4, FUN = "/")
  log_norm_mat <- log1p(norm_mat)

  sce_out <- sce_filt[common_genes, ]
  assay(sce_out, "counts")   <- corrected
  assay(sce_out, "norm")     <- norm_mat
  assay(sce_out, "log_norm") <- log_norm_mat
  assay(sce_out, "X")        <- log_norm_mat

  message("soupx.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce_out, output_path)
  message("soupx.R: done")
}

main <- function() {
  opts <- parse_args_soupx()
  run_soupx(opts$filtered_h5ad, opts$raw_dir, opts$output)
}

if (!interactive()) main()
