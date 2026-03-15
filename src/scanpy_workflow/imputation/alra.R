#!/usr/bin/env Rscript
# alra.R — ALRA imputation wrapper for Nextflow DSL2
#
# Usage:
#   Rscript alra.R --input <adata.h5ad> --output <imputed.h5ad>
#
# Input:  .h5ad with assay "log_norm" (log1p-normalized, genes x cells)
# Output: .h5ad with assay "alra_imputed" (genes x cells) and "X" updated.
#         Assays "counts", "norm", "log_norm" preserved unchanged.
#
# Requires ALRA from GitHub:
#   Rscript -e "remotes::install_github('nalab-stanford/ALRA')"
#   (see: make install-r-github)
#
# Note: ALRA expects cells x genes; this wrapper transposes before/after.

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(ALRA)
  library(Matrix)
})

parse_args_alra <- function() {
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

run_alra <- function(input_path, output_path) {
  message("alra.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!"log_norm" %in% assayNames(sce)) {
    stop("Input .h5ad must contain assay 'log_norm'")
  }

  # assay() -> genes x cells; ALRA expects cells x genes
  log_norm_cg <- t(as.matrix(assay(sce, "log_norm")))  # cells x genes

  message("alra.R: running ALRA imputation (",
          nrow(log_norm_cg), " cells x ", ncol(log_norm_cg), " genes)")
  alra_result <- alra(log_norm_cg)
  # alra() returns a 3-element list:
  #   [[1]] A_norm_rank_k      -- raw rank-k SVD approx (may contain negatives)
  #   [[2]] A_norm_rank_k_cor  -- threshold-corrected imputed matrix (non-negative)
  #   [[3]] d                  -- singular values
  # Use [[2]] for the biologically meaningful zero-restored output.
  imputed_cg  <- alra_result[[2]]   # A_norm_rank_k_cor

  # Transpose back to genes x cells
  imputed_gc <- t(imputed_cg)
  rownames(imputed_gc) <- rownames(sce)
  colnames(imputed_gc) <- colnames(sce)

  assay(sce, "alra_imputed") <- imputed_gc
  assay(sce, "X")            <- imputed_gc  # update adata.X

  message("alra.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("alra.R: done")
}

main <- function() {
  opts <- parse_args_alra()
  run_alra(opts$input, opts$output)
}

if (!interactive()) main()
