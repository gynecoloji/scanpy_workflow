#!/usr/bin/env Rscript
# tests/unit/r/test_mast.R — Unit tests for mast.R DE wrapper

suppressPackageStartupMessages({
  library(testthat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

# ── helpers ──────────────────────────────────────────────────────────────────

mast_r <- file.path(
  Sys.getenv("MAST_R_PATH",
             normalizePath(file.path(dirname(sys.frame(1)$ofile), "..",
                                     "..", "..", "src", "scanpy_workflow",
                                     "de", "mast.R"),
                           mustWork = FALSE)),
  fsep = "/"
)
if (!file.exists(mast_r)) {
  # fallback: relative to repo root
  mast_r <- normalizePath("src/scanpy_workflow/de/mast.R", mustWork = FALSE)
}
source(mast_r, local = TRUE)

# Build a small synthetic SCE with log_norm assay and leiden labels
make_fixture_sce <- function(n_genes = 50L, n_cells = 60L, n_clusters = 3L,
                              seed = 42L) {
  set.seed(seed)
  mat <- matrix(
    rpois(n_genes * n_cells, lambda = 1),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("Gene", seq_len(n_genes)),
                    paste0("Cell", seq_len(n_cells)))
  )
  mat_log <- log1p(mat / colSums(mat) * 1e4)
  cluster_labels <- rep(paste0("c", seq_len(n_clusters)),
                        length.out = n_cells)
  sce <- SingleCellExperiment(
    assays  = list(log_norm = mat_log),
    colData = data.frame(leiden = cluster_labels,
                         row.names = colnames(mat))
  )
  sce
}

write_fixture_h5ad <- function(sce, tmp_dir) {
  h5ad_path <- file.path(tmp_dir, "test_input.h5ad")
  writeH5AD(sce, h5ad_path)
  h5ad_path
}

# ── tests ─────────────────────────────────────────────────────────────────────

test_that("run_mast produces CSV with required columns", {
  skip_if_not_installed("MAST")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_out.csv")

  run_mast(h5ad, groupby = "leiden", output_path = out)

  expect_true(file.exists(out))
  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_true(all(c("group", "gene", "score", "logfoldchange", "pval", "pval_adj") %in%
                    names(df)))
})

test_that("run_mast group values match leiden labels", {
  skip_if_not_installed("MAST")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_groups.csv")

  run_mast(h5ad, groupby = "leiden", output_path = out)

  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_setequal(unique(df$group), c("c1", "c2", "c3"))
})

test_that("run_mast pval is in [0, 1]", {
  skip_if_not_installed("MAST")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_pval.csv")

  run_mast(h5ad, groupby = "leiden", output_path = out)

  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_true(all(df$pval     >= 0 & df$pval     <= 1))
  expect_true(all(df$pval_adj >= 0 & df$pval_adj <= 1))
})

test_that("run_mast pval_adj >= pval (BH is non-decreasing)", {
  skip_if_not_installed("MAST")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_bh.csv")

  run_mast(h5ad, groupby = "leiden", output_path = out)

  df <- read.csv(out, stringsAsFactors = FALSE)
  # BH adjusted p-value is always >= raw p-value (for non-trivial cases)
  expect_true(all(df$pval_adj + 1e-9 >= df$pval))
})

test_that("run_mast stops when log_norm assay is missing", {
  skip_if_not_installed("MAST")
  tmp <- tempdir()
  sce <- make_fixture_sce()
  # Replace log_norm with different name
  assayNames(sce) <- "wrong_assay"
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_fail.csv")

  expect_error(run_mast(h5ad, groupby = "leiden", output_path = out),
               regexp = "log_norm")
})

test_that("run_mast stops when groupby column is missing", {
  skip_if_not_installed("MAST")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "mast_fail2.csv")

  expect_error(run_mast(h5ad, groupby = "nonexistent_col", output_path = out),
               regexp = "colData missing")
})

# ── run ───────────────────────────────────────────────────────────────────────
if (!interactive()) {
  test_results <- testthat::test_file(
    normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
    reporter = "progress"
  )
}
