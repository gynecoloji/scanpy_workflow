# tests/unit/r/helper_fixture.R
# Shared helpers for R unit tests.
# Source this file at the top of each test file:
#   source(testthat::test_path("helper_fixture.R"))

library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)

#' Absolute path to the repository root, resolved via testthat::test_path().
#' Robust whether called from test_file(), test_dir(), or interactively.
project_root <- function() {
  normalizePath(file.path(testthat::test_path("."), "../../.."))
}

#' Absolute path to an R script under src/scanpy_workflow/.
#' @param ... Path components relative to src/scanpy_workflow/
script_path <- function(...) {
  file.path(project_root(), "src", "scanpy_workflow", ...)
}

#' Create a minimal synthetic .h5ad fixture.
#'
#' Writes an AnnData-compatible .h5ad with:
#'   - assay "X"        = log1p-normalized counts
#'   - assay "counts"   = raw integer counts
#'   - assay "norm"     = library-size normalized (pre-log)
#'   - assay "log_norm" = log1p(norm)
#'   - colData column "leiden" = random cluster assignments ("0" or "1")
#'
#' @param path    Output .h5ad path (character).
#' @param n_cells Integer; number of cells. Default 60.
#' @param n_genes Integer; number of genes. Default 80.
#' @return Invisibly, the SingleCellExperiment written.
make_fixture_h5ad <- function(path, n_cells = 60L, n_genes = 80L) {
  set.seed(42)
  counts_mat <- matrix(
    rpois(n_cells * n_genes, lambda = 3),
    nrow = n_genes, ncol = n_cells
  )
  # Ensure no all-zero cells (would break size-factor estimation in scran)
  counts_mat[1, counts_mat[1, ] == 0] <- 1L
  storage.mode(counts_mat) <- "integer"
  rownames(counts_mat) <- paste0("gene_", seq_len(n_genes))
  colnames(counts_mat) <- paste0("cell_", seq_len(n_cells))

  lib_sizes    <- colSums(counts_mat)
  norm_mat     <- sweep(counts_mat, 2, lib_sizes / 1e4, FUN = "/")
  log_norm_mat <- log1p(norm_mat)

  sce <- SingleCellExperiment(
    assays = list(
      X        = log_norm_mat,
      counts   = counts_mat,
      norm     = norm_mat,
      log_norm = log_norm_mat
    ),
    colData = data.frame(
      leiden = sample(c("0", "1"), n_cells, replace = TRUE),
      row.names = colnames(counts_mat)
    )
  )
  writeH5AD(sce, path)
  invisible(sce)
}

#' Create a minimal synthetic raw CellRanger-like matrix directory.
#'
#' Writes matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz into `dir`.
#' The raw dir has more barcodes (empty droplets) than the filtered data.
#' Gene names match those in make_fixture_h5ad() for the same n_genes.
#' Uses base R gzip connections — no external packages required.
#'
#' @param dir         Directory to create and populate.
#' @param n_cells_raw Number of barcodes in raw (including empty droplets). Default 300.
#' @param n_genes     Number of genes. Default 80. Must match fixture.
make_raw_cellranger_dir <- function(dir, n_cells_raw = 300L, n_genes = 80L) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(99)
  # Raw matrix: mostly sparse (empty droplets)
  counts_raw <- matrix(
    rpois(n_genes * n_cells_raw, lambda = 0.05),
    nrow = n_genes, ncol = n_cells_raw
  )
  # First 60 barcodes are "real" cells with higher counts
  counts_raw[, 1:60] <- counts_raw[, 1:60] + matrix(
    rpois(n_genes * 60, lambda = 3), nrow = n_genes
  )
  storage.mode(counts_raw) <- "integer"

  barcodes_raw <- paste0("cell_", seq_len(n_cells_raw), "-1")
  features     <- paste0("gene_", seq_len(n_genes))
  # Barcodes in filtered h5ad are "cell_1" through "cell_60"
  # SoupX requires the raw matrix to contain all filtered barcodes.
  # Rename first 60 raw barcodes to match filtered barcodes exactly.
  barcodes_raw[1:60] <- paste0("cell_", 1:60)

  mtx    <- Matrix::Matrix(counts_raw, sparse = TRUE)

  # matrix.mtx.gz — use base R gzip connection (no R.utils needed)
  con_gz <- gzfile(file.path(dir, "matrix.mtx.gz"), "w")
  Matrix::writeMM(mtx, con_gz)
  close(con_gz)

  # barcodes.tsv.gz
  con_bc <- gzfile(file.path(dir, "barcodes.tsv.gz"), "w")
  writeLines(barcodes_raw, con_bc)
  close(con_bc)

  # features.tsv.gz  (id <tab> name <tab> feature_type)
  con_ft <- gzfile(file.path(dir, "features.tsv.gz"), "w")
  writeLines(paste(features, features, "Gene Expression", sep = "\t"), con_ft)
  close(con_ft)

  invisible(dir)
}
