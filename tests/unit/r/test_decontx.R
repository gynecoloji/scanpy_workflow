# tests/unit/r/test_decontx.R
library(testthat)
library(zellkonverter)
library(SingleCellExperiment)

source(testthat::test_path("helper_fixture.R"))

run_decontx <- function(tmp_in, tmp_out) {
  ret <- system2(
    "Rscript",
    args = c(script_path("ambient", "decontx.R"),
             "--input",  tmp_in,
             "--output", tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_equal(attr(ret, "status"), 0,
    info = paste("decontx.R exited non-zero:\n", paste(ret, collapse = "\n")))
}

test_that("decontx.R produces output .h5ad with corrected counts and required layers", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in, n_cells = 60L, n_genes = 80L)
  run_decontx(tmp_in, tmp_out)

  sce_out <- readH5AD(tmp_out)
  expect_true("counts"   %in% assayNames(sce_out))
  expect_true("norm"     %in% assayNames(sce_out))
  expect_true("log_norm" %in% assayNames(sce_out))
  expect_true("X"        %in% assayNames(sce_out))

  counts_out <- as.matrix(assay(sce_out, "counts"))
  expect_true(all(counts_out >= 0), "corrected counts must be non-negative")

  X_mat        <- as.matrix(assay(sce_out, "X"))
  log_norm_mat <- as.matrix(assay(sce_out, "log_norm"))
  expect_true(isTRUE(all.equal(X_mat, log_norm_mat, tolerance = 1e-5)),
    info = "X must equal log_norm")

  sce_in <- readH5AD(tmp_in)
  expect_equal(nrow(sce_out), nrow(sce_in))
  expect_equal(ncol(sce_out), ncol(sce_in))
})

test_that("decontx.R adds decontX_contamination column to colData", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in)
  run_decontx(tmp_in, tmp_out)

  sce_out <- readH5AD(tmp_out)
  expect_true("decontX_contamination" %in% names(colData(sce_out)),
    info = "colData must contain decontX_contamination column")
})
