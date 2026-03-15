# tests/unit/r/test_soupx.R
library(testthat)
library(zellkonverter)
library(SingleCellExperiment)

source(testthat::test_path("helper_fixture.R"))

run_soupx <- function(tmp_in, tmp_raw, tmp_out) {
  ret <- system2(
    "Rscript",
    args = c(script_path("ambient", "soupx.R"),
             "--filtered-h5ad", tmp_in,
             "--raw-dir",       tmp_raw,
             "--output",        tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_equal(attr(ret, "status"), 0,
    info = paste("soupx.R exited non-zero:\n", paste(ret, collapse = "\n")))
}

test_that("soupx.R produces output .h5ad with corrected counts and required layers", {
  tmp_in  <- tempfile(fileext = ".h5ad")
  tmp_raw <- tempfile()
  tmp_out <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_raw, recursive = TRUE); unlink(tmp_out) },
          add = TRUE)

  make_fixture_h5ad(tmp_in, n_cells = 60L, n_genes = 80L)
  make_raw_cellranger_dir(tmp_raw, n_cells_raw = 300L, n_genes = 80L)
  run_soupx(tmp_in, tmp_raw, tmp_out)

  sce_out <- readH5AD(tmp_out)
  expect_true("counts"   %in% assayNames(sce_out))
  expect_true("norm"     %in% assayNames(sce_out))
  expect_true("log_norm" %in% assayNames(sce_out))
  expect_true("X"        %in% assayNames(sce_out))

  counts_out <- as.matrix(assay(sce_out, "counts"))
  expect_true(all(counts_out >= 0))

  X_mat        <- as.matrix(assay(sce_out, "X"))
  log_norm_mat <- as.matrix(assay(sce_out, "log_norm"))
  expect_true(isTRUE(all.equal(X_mat, log_norm_mat, tolerance = 1e-5)),
    info = "X must equal log_norm")

  # Cell and gene count preserved
  sce_in <- readH5AD(tmp_in)
  expect_equal(ncol(sce_out), ncol(sce_in))
  expect_equal(nrow(sce_out), nrow(sce_in))
})

test_that("soupx.R exits non-zero when raw dir is missing filtered barcodes", {
  tmp_in       <- tempfile(fileext = ".h5ad")
  tmp_raw_bad  <- tempfile()
  tmp_out      <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_raw_bad, recursive = TRUE); unlink(tmp_out) },
          add = TRUE)

  # Fixture: 60 filtered cells named "cell_1" through "cell_60"
  make_fixture_h5ad(tmp_in, n_cells = 60L, n_genes = 80L)
  # Raw dir with completely different barcode names (no overlap with filtered)
  make_raw_cellranger_dir(tmp_raw_bad, n_cells_raw = 300L, n_genes = 80L)
  # Overwrite barcodes.tsv.gz with non-matching names
  con_bc <- gzfile(file.path(tmp_raw_bad, "barcodes.tsv.gz"), "w")
  writeLines(paste0("unrelated_bc_", seq_len(300), "-1"), con_bc)
  close(con_bc)

  ret <- system2(
    "Rscript",
    args = c(script_path("ambient", "soupx.R"),
             "--filtered-h5ad", tmp_in,
             "--raw-dir",       tmp_raw_bad,
             "--output",        tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_true(attr(ret, "status") != 0,
    info = "soupx.R must exit non-zero when filtered barcodes are absent from raw dir")
})
