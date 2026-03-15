#!/usr/bin/env Rscript
# decontx.R — DecontX ambient RNA removal wrapper for Nextflow DSL2
#
# Usage:
#   Rscript decontx.R --input <filtered.h5ad> --output <decontaminated.h5ad>
#
# Uses bioconductor-celda::decontX. Passes leiden cluster labels (z) when
# present in colData for improved contamination estimation.
# Writes corrected integer counts (rounded) + recomputed norm/log_norm/X.
# Adds colData column "decontX_contamination".

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(celda)
  library(Matrix)
})

parse_args_decontx <- function() {
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

run_decontx <- function(input_path, output_path) {
  message("decontx.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!"counts" %in% assayNames(sce)) {
    stop("Input .h5ad must contain assay 'counts'")
  }
  counts_mat <- assay(sce, "counts")

  # Pass cluster labels when available — improves contamination estimation
  z_arg <- NULL
  if ("leiden" %in% names(colData(sce))) {
    z_arg <- as.integer(factor(colData(sce)$leiden))
    message("decontx.R: using leiden cluster labels (", length(unique(z_arg)), " clusters)")
  } else {
    message("decontx.R: no cluster labels found; decontX will estimate internally")
  }

  message("decontx.R: running decontX")
  sce_temp    <- SingleCellExperiment(assays = list(counts = counts_mat))
  decontx_res <- if (!is.null(z_arg)) decontX(sce_temp, z = z_arg) else decontX(sce_temp)

  # Corrected counts (rounded to integers, floored at 0)
  corrected <- round(decontXcounts(decontx_res))
  corrected[corrected < 0L] <- 0L
  storage.mode(corrected) <- "integer"

  # Recompute downstream layers from corrected counts
  lib_sizes    <- colSums(corrected)
  lib_sizes[lib_sizes == 0] <- 1
  norm_mat     <- sweep(corrected, 2, lib_sizes / 1e4, FUN = "/")
  log_norm_mat <- log1p(norm_mat)

  assay(sce, "counts")   <- corrected
  assay(sce, "norm")     <- norm_mat
  assay(sce, "log_norm") <- log_norm_mat
  assay(sce, "X")        <- log_norm_mat

  sce$decontX_contamination <- colData(decontx_res)$decontX_contamination

  message("decontx.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("decontx.R: done")
}

main <- function() {
  opts <- parse_args_decontx()
  run_decontx(opts$input, opts$output)
}

if (!interactive()) main()
