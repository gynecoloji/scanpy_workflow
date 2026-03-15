# tests/unit/r/test_sctype.R
library(testthat)
library(zellkonverter)
library(SingleCellExperiment)
library(jsonlite)

source(testthat::test_path("helper_fixture.R"))

make_markers_json <- function(path, gene_names) {
  markers <- list(
    `T cells` = list(pos = gene_names[1:3], neg = list()),
    `B cells` = list(pos = gene_names[4:5], neg = gene_names[1:2])
  )
  writeLines(toJSON(markers, auto_unbox = TRUE, pretty = TRUE), path)
}

run_sctype <- function(tmp_in, tmp_markers, tmp_out, groupby = "leiden") {
  ret <- system2(
    "Rscript",
    args = c(script_path("annotation", "sctype.R"),
             "--input",   tmp_in,
             "--markers", tmp_markers,
             "--groupby", groupby,
             "--output",  tmp_out),
    stdout = TRUE, stderr = TRUE
  )
  expect_equal(attr(ret, "status"), 0,
    info = paste("sctype.R exited non-zero:\n", paste(ret, collapse = "\n")))
}

test_that("sctype.R adds sctype_annotation column to colData", {
  tmp_in      <- tempfile(fileext = ".h5ad")
  tmp_markers <- tempfile(fileext = ".json")
  tmp_out     <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_markers); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in, n_cells = 60L, n_genes = 80L)
  make_markers_json(tmp_markers, paste0("gene_", 1:80))
  run_sctype(tmp_in, tmp_markers, tmp_out)

  sce_out     <- readH5AD(tmp_out)
  expect_true("sctype_annotation" %in% names(colData(sce_out)),
    info = "colData must contain 'sctype_annotation'")

  annotations <- colData(sce_out)$sctype_annotation
  expect_true(all(!is.na(annotations) & annotations != ""))

  valid_labels <- c("T cells", "B cells", "Unknown")
  expect_true(all(annotations %in% valid_labels),
    info = paste("Unexpected labels:", paste(setdiff(unique(annotations), valid_labels), collapse = ", ")))

  sce_in <- readH5AD(tmp_in)
  expect_equal(ncol(sce_out), ncol(sce_in))
  expect_equal(nrow(sce_out), nrow(sce_in))
})

test_that("sctype.R assigns same annotation to all cells in the same cluster", {
  tmp_in      <- tempfile(fileext = ".h5ad")
  tmp_markers <- tempfile(fileext = ".json")
  tmp_out     <- tempfile(fileext = ".h5ad")
  on.exit({ unlink(tmp_in); unlink(tmp_markers); unlink(tmp_out) }, add = TRUE)

  make_fixture_h5ad(tmp_in)
  make_markers_json(tmp_markers, paste0("gene_", 1:80))
  run_sctype(tmp_in, tmp_markers, tmp_out)

  sce_out     <- readH5AD(tmp_out)
  cd          <- as.data.frame(colData(sce_out))
  # All cells in the same cluster get the same annotation
  per_cluster <- tapply(cd$sctype_annotation, cd$leiden,
                        function(x) length(unique(x)))
  expect_true(all(per_cluster == 1),
    info = "All cells in the same cluster must have the same sctype_annotation")
})
