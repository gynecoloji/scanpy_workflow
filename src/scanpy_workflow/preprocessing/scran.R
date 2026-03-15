#!/usr/bin/env Rscript
# scran.R — scran normalization wrapper for Nextflow DSL2
#
# Usage:
#   Rscript scran.R --input <in.h5ad> --output <out.h5ad>
#
# Layer contract (identical to normalize.py library_size path):
#   assay "counts"   — raw integer counts (unchanged)
#   assay "norm"     — scran size-factor normalized, pre-log
#   assay "log_norm" — log1p(norm)
#   assay "X"        — set to log_norm (adata.X = layers["log_norm"])
#
# Requires: bioconductor-scran >= 1.30, bioconductor-scater >= 1.30

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(Matrix)
})

parse_args_scran <- function() {
  option_list <- list(
    make_option("--input",  type = "character", help = "Input .h5ad path"),
    make_option("--output", type = "character", help = "Output .h5ad path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$input) || is.null(opts$output)) {
    print_help(parser)
    stop("--input and --output are required", call. = FALSE)
  }
  opts
}

run_scran <- function(input_path, output_path) {
  message("scran.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!"counts" %in% assayNames(sce)) {
    stop("Input .h5ad must have assay 'counts' (raw integer counts)")
  }
  counts_mat <- assay(sce, "counts")

  message("scran.R: estimating size factors via pooling")
  # quickCluster → computeSumFactors for robust size-factor estimation
  clusters <- quickCluster(counts_mat, min.size = 10)
  sce_temp <- SingleCellExperiment(assays = list(counts = counts_mat))
  sce_temp <- computeSumFactors(sce_temp, cluster = clusters)
  # logNormCounts(log = FALSE) stores size-factor-normalized counts in "normcounts"
  # Requires scater >= 1.18; pinned to >= 1.30 in env_r.yaml
  sce_temp <- logNormCounts(sce_temp, log = FALSE)

  norm_mat     <- normcounts(sce_temp)
  log_norm_mat <- log1p(norm_mat)

  # Write all three layers + X back to the original SCE (preserving metadata)
  assay(sce, "counts")   <- counts_mat
  assay(sce, "norm")     <- norm_mat
  assay(sce, "log_norm") <- log_norm_mat
  assay(sce, "X")        <- log_norm_mat   # adata.X = layers["log_norm"]

  message("scran.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("scran.R: done")
}

main <- function() {
  opts <- parse_args_scran()
  run_scran(opts$input, opts$output)
}

if (!interactive()) main()
