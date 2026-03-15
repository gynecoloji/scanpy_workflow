# tests/unit/r/test_scran.R
library(testthat)
library(zellkonverter)
library(SingleCellExperiment)

source(testthat::test_path("helper_fixture.R"))

run_scran <- function(tmp_in, tmp_out) {
  ret <- system2(
    "Rscript",
    args = c(script_path("preprocessing", "scran.R"),
             "--input",  tmp_in,
             "--output", tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_equal(attr(ret, "status"), 0,
    info = paste("scran.R exited non-zero:\n", paste(ret, collapse = "\n")))
}

test_that("scran.R writes three required layers and sets X = log_norm", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in)
  run_scran(tmp_in, tmp_out)

  sce_out <- readH5AD(tmp_out)

  expect_true("counts"   %in% assayNames(sce_out), "Missing assay: counts")
  expect_true("norm"     %in% assayNames(sce_out), "Missing assay: norm")
  expect_true("log_norm" %in% assayNames(sce_out), "Missing assay: log_norm")
  expect_true("X"        %in% assayNames(sce_out), "Missing assay: X")

  # X must equal log_norm (downstream Python depends on this)
  X_mat        <- as.matrix(assay(sce_out, "X"))
  log_norm_mat <- as.matrix(assay(sce_out, "log_norm"))
  expect_true(isTRUE(all.equal(X_mat, log_norm_mat, tolerance = 1e-5)),
    info = "assay 'X' must equal assay 'log_norm'")

  # counts must be unchanged
  sce_in    <- readH5AD(tmp_in)
  expect_equal(as.matrix(assay(sce_in, "counts")),
               as.matrix(assay(sce_out, "counts")),
               info = "layers['counts'] must not change")

  # log_norm must be log1p(norm)
  norm_mat <- as.matrix(assay(sce_out, "norm"))
  expect_true(isTRUE(all.equal(log_norm_mat, log1p(norm_mat), tolerance = 1e-5)),
    info = "log_norm must be log1p(norm)")
})

test_that("scran.R preserves cell and gene names", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in)
  run_scran(tmp_in, tmp_out)

  sce_in  <- readH5AD(tmp_in)
  sce_out <- readH5AD(tmp_out)
  expect_equal(rownames(sce_in), rownames(sce_out))
  expect_equal(colnames(sce_in), colnames(sce_out))
})
