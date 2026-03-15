# tests/unit/r/test_alra.R
library(testthat)
library(zellkonverter)
library(SingleCellExperiment)

source(testthat::test_path("helper_fixture.R"))

run_alra <- function(tmp_in, tmp_out) {
  ret <- system2(
    "Rscript",
    args = c(script_path("imputation", "alra.R"),
             "--input",  tmp_in,
             "--output", tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_equal(attr(ret, "status"), 0,
    info = paste("alra.R exited non-zero:\n", paste(ret, collapse = "\n")))
}

test_that("alra.R imputes dropout and updates X assay", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)

  # 200 cells for reliable low-rank structure (avoids ALRA thresholding flakiness)
  make_fixture_h5ad(tmp_in, n_cells = 200L, n_genes = 80L)
  run_alra(tmp_in, tmp_out)

  sce_out <- readH5AD(tmp_out)

  expect_true("counts"       %in% assayNames(sce_out))
  expect_true("norm"         %in% assayNames(sce_out))
  expect_true("log_norm"     %in% assayNames(sce_out))
  expect_true("alra_imputed" %in% assayNames(sce_out))
  expect_true("X"            %in% assayNames(sce_out))

  sce_in <- readH5AD(tmp_in)
  expect_equal(nrow(sce_out), nrow(sce_in))
  expect_equal(ncol(sce_out), ncol(sce_in))

  X_out       <- as.matrix(assay(sce_out, "X"))
  log_norm_in <- as.matrix(assay(sce_in, "log_norm"))

  # Imputed X must be non-negative (ALRA [[2]] is threshold-corrected)
  expect_true(all(X_out >= 0),
    label = "imputed X must be non-negative (ALRA A_norm_rank_k_cor)")

  # ALRA should reduce the zero fraction in X relative to log_norm
  expect_lt(sum(X_out == 0), sum(log_norm_in == 0),
    label = "ALRA should reduce zero fraction in X")

  # Original counts must not be modified
  expect_equal(as.matrix(assay(sce_in, "counts")),
               as.matrix(assay(sce_out, "counts")))
})
