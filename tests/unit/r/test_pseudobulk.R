#!/usr/bin/env Rscript
# tests/unit/r/test_pseudobulk.R — Unit tests for pseudobulk.R DE wrapper

suppressPackageStartupMessages({
  library(testthat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})

# ── helpers ──────────────────────────────────────────────────────────────────

pseudobulk_r <- file.path(
  Sys.getenv("PSEUDOBULK_R_PATH",
             normalizePath(file.path(dirname(sys.frame(1)$ofile), "..",
                                     "..", "..", "src", "scanpy_workflow",
                                     "de", "pseudobulk.R"),
                           mustWork = FALSE)),
  fsep = "/"
)
if (!file.exists(pseudobulk_r)) {
  pseudobulk_r <- normalizePath("src/scanpy_workflow/de/pseudobulk.R",
                                 mustWork = FALSE)
}
source(pseudobulk_r, local = TRUE)

# Build a synthetic SCE suitable for pseudobulk:
#   - 2 samples × 40 cells each = 80 cells total
#   - 3 clusters assigned in blocks so each (cluster, sample) has ≥10 cells
make_fixture_sce <- function(n_genes = 30L, n_cells_per_sample = 40L,
                              n_clusters = 3L, seed = 42L) {
  set.seed(seed)
  n_cells  <- n_cells_per_sample * 2L
  counts   <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("Gene", seq_len(n_genes)),
                    paste0("Cell", seq_len(n_cells)))
  )
  # Cluster labels: first 2/3 cells per sample to cluster 1, rest split
  clusters <- rep(c(rep("c1", 14L), rep("c2", 13L), rep("c3", 13L)),
                  times = 2L)
  samples  <- rep(c("sampleA", "sampleB"), each = n_cells_per_sample)

  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    colData = data.frame(leiden = clusters,
                         sample = samples,
                         row.names = colnames(counts))
  )
  sce
}

write_fixture_h5ad <- function(sce, tmp_dir) {
  h5ad_path <- file.path(tmp_dir, "pb_test_input.h5ad")
  writeH5AD(sce, h5ad_path)
  h5ad_path
}

# ── tests ─────────────────────────────────────────────────────────────────────

test_that("run_pseudobulk produces CSV with required columns", {
  skip_if_not_installed("DESeq2")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_out.csv")

  run_pseudobulk(h5ad, groupby = "leiden", sample_key = "sample",
                 min_cells = 10L, output_path = out)

  expect_true(file.exists(out))
  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_true(all(c("group", "gene", "score", "logfoldchange",
                    "pval", "pval_adj") %in% names(df)))
})

test_that("run_pseudobulk group values match leiden labels", {
  skip_if_not_installed("DESeq2")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_groups.csv")

  run_pseudobulk(h5ad, groupby = "leiden", sample_key = "sample",
                 min_cells = 10L, output_path = out)

  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_setequal(unique(df$group), c("c1", "c2", "c3"))
})

test_that("run_pseudobulk pval is in [0, 1]", {
  skip_if_not_installed("DESeq2")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_pval.csv")

  run_pseudobulk(h5ad, groupby = "leiden", sample_key = "sample",
                 min_cells = 10L, output_path = out)

  df <- read.csv(out, stringsAsFactors = FALSE)
  expect_true(all(df$pval     >= 0 & df$pval     <= 1))
  expect_true(all(df$pval_adj >= 0 & df$pval_adj <= 1))
})

test_that("aggregate_pseudobulk sums counts correctly", {
  counts <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
                   dimnames = list(c("G1", "G2"), c("A", "B", "C")))
  clusters <- c("cl1", "cl1", "cl2")
  samples  <- c("s1",  "s1",  "s1")
  pb <- aggregate_pseudobulk(counts, clusters, samples)
  expect_equal(pb$mat["G1", "cl1__s1"], 1 + 3)   # cells A + B
  expect_equal(pb$mat["G2", "cl1__s1"], 2 + 4)
  expect_equal(pb$mat["G1", "cl2__s1"], 5)
})

test_that("run_pseudobulk stops when counts assay is missing", {
  skip_if_not_installed("DESeq2")
  tmp <- tempdir()
  sce <- make_fixture_sce()
  assayNames(sce) <- "wrong_assay"
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_fail.csv")

  expect_error(
    run_pseudobulk(h5ad, groupby = "leiden", sample_key = "sample",
                   min_cells = 10L, output_path = out),
    regexp = "counts"
  )
})

test_that("run_pseudobulk stops when groupby column is missing", {
  skip_if_not_installed("DESeq2")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_fail2.csv")

  expect_error(
    run_pseudobulk(h5ad, groupby = "bad_col", sample_key = "sample",
                   min_cells = 10L, output_path = out),
    regexp = "colData missing"
  )
})

test_that("run_pseudobulk stops when sample_key column is missing", {
  skip_if_not_installed("DESeq2")
  tmp  <- tempdir()
  sce  <- make_fixture_sce()
  h5ad <- write_fixture_h5ad(sce, tmp)
  out  <- file.path(tmp, "pb_fail3.csv")

  expect_error(
    run_pseudobulk(h5ad, groupby = "leiden", sample_key = "bad_sample",
                   min_cells = 10L, output_path = out),
    regexp = "colData missing"
  )
})

# ── run ───────────────────────────────────────────────────────────────────────
if (!interactive()) {
  test_results <- testthat::test_file(
    normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
    reporter = "progress"
  )
}
