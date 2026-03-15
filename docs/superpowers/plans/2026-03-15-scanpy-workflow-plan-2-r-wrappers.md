# Scanpy Workflow — Plan 2: R Wrappers

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement all five R wrapper scripts (`soupx.R`, `decontx.R`, `scran.R`, `alra.R`, `sctype.R`) with `testthat` unit tests, the `env_r` conda environment YAML, and CI wiring.

**Architecture:** Each R script is a standalone CLI (`Rscript script.R --input ... --output ...`) called directly by Nextflow. Each script defines a `main()` function and guards it with `if (!interactive()) main()` so tests can `source()` and call functions directly. All `.h5ad` I/O uses `zellkonverter`. Argument parsing uses `optparse`. All tests use `testthat::test_path()` for path resolution (robust under both `test_file()` and `test_dir()` invocations).

**Tech Stack:** R ≥ 4.3, zellkonverter, SingleCellExperiment, optparse, SoupX, bioconductor-celda (DecontX), scran (≥1.30), scater (≥1.30), ALRA (GitHub), scType (GitHub), testthat, Matrix

**Spec:** `docs/superpowers/specs/2026-03-15-scanpy-workflow-design.md`

---

## Chunk 1: Setup — env_r.yaml + test helpers

### Task 1: env_r conda environment YAML

**Files:**
- Create: `envs/env_r.yaml`

- [ ] **Step 1: Write `envs/env_r.yaml`**

```yaml
# envs/env_r.yaml
name: env_r
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - r-base=4.3
  - r-optparse
  - r-matrix
  - r-testthat
  - r-jsonlite
  - bioconductor-singlecellexperiment
  - bioconductor-zellkonverter
  - bioconductor-scran>=1.30        # logNormCounts(log=FALSE) requires >=1.18; pin >=1.30 (Bioc 3.18+)
  - bioconductor-scater>=1.30       # normcounts() accessor
  - bioconductor-celda              # provides decontX
  - r-soupx
  # ALRA and scType are not on conda; install from GitHub in post-install step:
  #   Rscript -e "remotes::install_github('nalab-stanford/ALRA')"
  #   Rscript -e "remotes::install_github('IanevskiAleksandr/sc-type')"
  # Run: make install-r-github  (see Makefile target below)
  - r-remotes                       # for GitHub installs above
  - r-biocmanager
```

- [ ] **Step 2: Add Makefile target for post-install GitHub packages**

Append to the project root `Makefile` (create if absent):

```makefile
## Install GitHub-only R packages (ALRA, scType) into the active R library.
## Run once after: conda env create -f envs/env_r.yaml && conda activate env_r
install-r-github:
	Rscript -e "remotes::install_github('nalab-stanford/ALRA',              upgrade = 'never')"
	Rscript -e "remotes::install_github('IanevskiAleksandr/sc-type',        upgrade = 'never')"
```

- [ ] **Step 3: Commit**

```bash
git add envs/env_r.yaml Makefile
git commit -m "feat: add env_r conda environment YAML and Makefile install target"
```

---

### Task 2: Shared test helper

**Files:**
- Create: `tests/unit/r/helper_fixture.R`

- [ ] **Step 1: Write `tests/unit/r/helper_fixture.R`**

This helper is `source()`d by all R test files. It creates a minimal synthetic `.h5ad` fixture on disk and provides shared path-resolution utilities.

```R
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
```

- [ ] **Step 2: Commit**

```bash
git add tests/unit/r/helper_fixture.R
git commit -m "feat: add R test fixture helper"
```

---

## Chunk 2: scran normalization

### Task 3: scran.R — normalization wrapper

**Files:**
- Create: `src/scanpy_workflow/preprocessing/scran.R`
- Create: `tests/unit/r/test_scran.R`

**Layer contract (must match `normalize.py`):**
- `assay(sce, "counts")` — raw integer counts (unchanged)
- `assay(sce, "norm")`   — scran size-factor normalized, pre-log
- `assay(sce, "log_norm")` — log1p of norm
- `assay(sce, "X")` — must equal `log_norm`

- [ ] **Step 1: Write the failing test**

```R
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_scran.R')"
```

Expected: FAIL — `src/scanpy_workflow/preprocessing/scran.R` does not exist.

- [ ] **Step 3: Implement `src/scanpy_workflow/preprocessing/scran.R`**

```R
#!/usr/bin/env Rscript
# scran.R — scran normalization wrapper for Nextflow DSL2
#
# Usage:
#   Rscript scran.R --input <in.h5ad> --output <out.h5ad>
#
# Layer contract (identical to normalize.py library_size path):
#   assay "counts"   — raw integer counts (unchanged)
#   assay "norm"     — scran size-factor normalized, pre-log
#   assay "log_norm" — log1p(norm)
#   assay "X"        — set to log_norm (adata.X = layers["log_norm"])
#
# Requires: bioconductor-scran >= 1.30, bioconductor-scater >= 1.30

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(Matrix)
})

parse_args_scran <- function() {
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

run_scran <- function(input_path, output_path) {
  message("scran.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!"counts" %in% assayNames(sce)) {
    stop("Input .h5ad must have assay 'counts' (raw integer counts)")
  }
  counts_mat <- assay(sce, "counts")

  message("scran.R: estimating size factors via pooling")
  # quickCluster → computeSumFactors for robust size-factor estimation
  clusters <- quickCluster(counts_mat, min.size = 10)
  sce_temp <- SingleCellExperiment(assays = list(counts = counts_mat))
  sce_temp <- computeSumFactors(sce_temp, cluster = clusters)
  # logNormCounts(log = FALSE) stores size-factor-normalized counts in "normcounts"
  # Requires scater >= 1.18; pinned to >= 1.30 in env_r.yaml
  sce_temp <- logNormCounts(sce_temp, log = FALSE)

  norm_mat     <- normcounts(sce_temp)
  log_norm_mat <- log1p(norm_mat)

  # Write all three layers + X back to the original SCE (preserving metadata)
  assay(sce, "counts")   <- counts_mat
  assay(sce, "norm")     <- norm_mat
  assay(sce, "log_norm") <- log_norm_mat
  assay(sce, "X")        <- log_norm_mat   # adata.X = layers["log_norm"]

  message("scran.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("scran.R: done")
}

main <- function() {
  opts <- parse_args_scran()
  run_scran(opts$input, opts$output)
}

if (!interactive()) main()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_scran.R')"
```

Expected: All 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/preprocessing/scran.R tests/unit/r/test_scran.R
git commit -m "feat: add scran normalization R wrapper with testthat tests"
```

---

## Chunk 3: Ambient RNA removal

### Task 4: soupx.R — SoupX ambient RNA removal

**Files:**
- Create: `src/scanpy_workflow/ambient/soupx.R`
- Create: `tests/unit/r/test_soupx.R`

**Interface:**
```
Rscript soupx.R \
  --filtered-h5ad <filtered.h5ad> \
  --raw-dir <cellranger_raw_feature_bc_matrix/> \
  --output <decontaminated.h5ad>
```

**Behavior:** Reads the raw CellRanger directory using a base R helper (no Seurat dependency). Corrected counts replace `assay(sce, "counts")`. Recomputes `norm`, `log_norm`, `X` from corrected counts.

- [ ] **Step 1: Write the failing test**

```R
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_soupx.R')"
```

Expected: FAIL — `src/scanpy_workflow/ambient/soupx.R` does not exist.

- [ ] **Step 3: Implement `src/scanpy_workflow/ambient/soupx.R`**

```R
#!/usr/bin/env Rscript
# soupx.R — SoupX ambient RNA removal wrapper for Nextflow DSL2
#
# Usage:
#   Rscript soupx.R \
#     --filtered-h5ad <filtered.h5ad> \
#     --raw-dir       <cellranger_raw_feature_bc_matrix/> \
#     --output        <decontaminated.h5ad>
#
# Reads the raw CellRanger directory with a base R helper (no Seurat required).
# Writes corrected integer counts + recomputed norm/log_norm/X layers.

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SoupX)
  library(Matrix)
})

parse_args_soupx <- function() {
  option_list <- list(
    make_option("--filtered-h5ad", type = "character", dest = "filtered_h5ad",
                help = "Filtered .h5ad input path"),
    make_option("--raw-dir", type = "character", dest = "raw_dir",
                help = "CellRanger raw_feature_bc_matrix directory"),
    make_option("--output", type = "character", help = "Output .h5ad path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$filtered_h5ad) || is.null(opts$raw_dir) || is.null(opts$output)) {
    print_help(parser)
    stop("--filtered-h5ad, --raw-dir, and --output are required", call. = FALSE)
  }
  opts
}

#' Read a 10x CellRanger matrix directory (gzipped files) without Seurat.
#' Returns a genes × barcodes sparse matrix (dgCMatrix).
read_10x_matrix <- function(dir) {
  bc_con   <- gzfile(file.path(dir, "barcodes.tsv.gz"), "r")
  barcodes <- readLines(bc_con)
  close(bc_con)

  ft_con     <- gzfile(file.path(dir, "features.tsv.gz"), "r")
  feat_lines <- readLines(ft_con)
  close(ft_con)
  feature_ids <- sub("\t.*", "", feat_lines)

  mtx_con <- gzfile(file.path(dir, "matrix.mtx.gz"), "r")
  mtx     <- Matrix::readMM(mtx_con)
  close(mtx_con)

  rownames(mtx) <- feature_ids
  colnames(mtx) <- barcodes
  as(mtx, "dgCMatrix")
}

run_soupx <- function(filtered_h5ad, raw_dir, output_path) {
  message("soupx.R: reading filtered .h5ad from ", filtered_h5ad)
  sce_filt <- readH5AD(filtered_h5ad, reader = "R")

  if (!"counts" %in% assayNames(sce_filt)) {
    stop("Filtered .h5ad must contain assay 'counts'")
  }
  toc <- assay(sce_filt, "counts")   # table of counts (filtered cells)

  message("soupx.R: reading raw CellRanger directory from ", raw_dir)
  tod <- read_10x_matrix(raw_dir)    # table of droplets (all barcodes incl. empty)

  # Align features: raw dir may have a different feature set
  common_genes <- intersect(rownames(tod), rownames(toc))
  if (length(common_genes) == 0) {
    stop("No common genes between raw dir and filtered .h5ad. Check gene name format.")
  }
  tod <- tod[common_genes, ]
  toc <- toc[common_genes, ]

  # Verify that raw matrix contains all filtered barcodes
  missing_bc <- setdiff(colnames(toc), colnames(tod))
  if (length(missing_bc) > 0) {
    stop(sprintf(
      "soupx.R: %d filtered barcodes not found in raw dir (e.g. %s). Check that --raw-dir points to the correct CellRanger raw_feature_bc_matrix.",
      length(missing_bc), paste(head(missing_bc, 3), collapse = ", ")
    ))
  }

  message("soupx.R: building SoupChannel and estimating contamination")
  sc <- SoupChannel(tod, toc)

  # Use leiden clusters if available; otherwise quick-cluster from lib sizes
  if ("leiden" %in% names(colData(sce_filt))) {
    cl <- setNames(as.character(colData(sce_filt)[colnames(toc), "leiden"]),
                   colnames(toc))
  } else {
    lib <- colSums(toc)
    cl  <- setNames(ifelse(lib > median(lib), "1", "0"), colnames(toc))
  }
  sc <- setClusters(sc, cl)
  sc <- autoEstCont(sc, doPlot = FALSE)

  message("soupx.R: adjusting counts")
  corrected <- adjustCounts(sc, roundToInt = TRUE)
  # adjustCounts may change column order; re-align to toc column order
  corrected <- corrected[, colnames(toc)]

  # Recompute downstream layers from corrected counts
  lib_sizes    <- colSums(corrected)
  lib_sizes[lib_sizes == 0] <- 1
  norm_mat     <- sweep(corrected, 2, lib_sizes / 1e4, FUN = "/")
  log_norm_mat <- log1p(norm_mat)

  sce_out <- sce_filt[common_genes, ]
  assay(sce_out, "counts")   <- corrected
  assay(sce_out, "norm")     <- norm_mat
  assay(sce_out, "log_norm") <- log_norm_mat
  assay(sce_out, "X")        <- log_norm_mat

  message("soupx.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce_out, output_path)
  message("soupx.R: done")
}

main <- function() {
  opts <- parse_args_soupx()
  run_soupx(opts$filtered_h5ad, opts$raw_dir, opts$output)
}

if (!interactive()) main()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_soupx.R')"
```

Expected: 1 test PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/ambient/soupx.R tests/unit/r/test_soupx.R
git commit -m "feat: add SoupX ambient RNA removal R wrapper with testthat tests"
```

---

### Task 5: decontx.R — DecontX ambient RNA removal

**Files:**
- Create: `src/scanpy_workflow/ambient/decontx.R`
- Create: `tests/unit/r/test_decontx.R`

**Interface:**
```
Rscript decontx.R --input <filtered.h5ad> --output <decontaminated.h5ad>
```

**Behavior:** DecontX (from `bioconductor-celda`) needs only filtered counts. Passes leiden cluster assignments (`z`) when available for improved estimation. Corrected counts replace `assay(sce, "counts")`. Recomputes `norm`, `log_norm`, `X`. Adds `colData$decontX_contamination`.

- [ ] **Step 1: Write the failing test**

```R
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_decontx.R')"
```

Expected: FAIL — `src/scanpy_workflow/ambient/decontx.R` does not exist.

- [ ] **Step 3: Implement `src/scanpy_workflow/ambient/decontx.R`**

```R
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_decontx.R')"
```

Expected: All 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/ambient/decontx.R tests/unit/r/test_decontx.R
git commit -m "feat: add DecontX ambient RNA removal R wrapper with testthat tests"
```

---

## Chunk 4: Imputation

### Task 6: alra.R — ALRA imputation wrapper

**Files:**
- Create: `src/scanpy_workflow/imputation/alra.R`
- Create: `tests/unit/r/test_alra.R`

**Interface:**
```
Rscript alra.R --input <adata.h5ad> --output <imputed.h5ad>
```

**Behavior:** ALRA imputes from log-normalized data (`assay(sce, "log_norm")`). The imputed matrix is stored as `assay(sce, "alra_imputed")` and also replaces `assay(sce, "X")`. Original `counts`, `norm`, `log_norm` are preserved.

> **Note:** ALRA expects cells × genes. `assay()` is genes × cells. Transpose before and after.

- [ ] **Step 1: Write the failing test**

```R
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_alra.R')"
```

Expected: FAIL — `src/scanpy_workflow/imputation/alra.R` does not exist.

- [ ] **Step 3: Implement `src/scanpy_workflow/imputation/alra.R`**

```R
#!/usr/bin/env Rscript
# alra.R — ALRA imputation wrapper for Nextflow DSL2
#
# Usage:
#   Rscript alra.R --input <adata.h5ad> --output <imputed.h5ad>
#
# Input:  .h5ad with assay "log_norm" (log1p-normalized, genes × cells)
# Output: .h5ad with assay "alra_imputed" (genes × cells) and "X" updated.
#         Assays "counts", "norm", "log_norm" preserved unchanged.
#
# Requires ALRA from GitHub:
#   Rscript -e "remotes::install_github('nalab-stanford/ALRA')"
#   (see: make install-r-github)
#
# Note: ALRA expects cells × genes; this wrapper transposes before/after.

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

  # assay() → genes × cells; ALRA expects cells × genes
  log_norm_cg <- t(as.matrix(assay(sce, "log_norm")))  # cells × genes

  message("alra.R: running ALRA imputation (",
          nrow(log_norm_cg), " cells × ", ncol(log_norm_cg), " genes)")
  alra_result <- alra(log_norm_cg)
  # alra() returns a 3-element list:
  #   [[1]] A_norm_rank_k      — raw rank-k SVD approx (may contain negatives)
  #   [[2]] A_norm_rank_k_cor  — threshold-corrected imputed matrix (non-negative)
  #   [[3]] d                  — singular values
  # Use [[2]] for the biologically meaningful zero-restored output.
  imputed_cg  <- alra_result[[2]]   # A_norm_rank_k_cor

  # Transpose back to genes × cells
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_alra.R')"
```

Expected: 1 test PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/imputation/alra.R tests/unit/r/test_alra.R
git commit -m "feat: add ALRA imputation R wrapper with testthat tests"
```

---

## Chunk 5: Annotation

### Task 7: sctype.R — scType automated annotation wrapper

**Files:**
- Create: `src/scanpy_workflow/annotation/sctype.R`
- Create: `tests/unit/r/test_sctype.R`

**Interface:**
```
Rscript sctype.R \
  --input   <clustered.h5ad> \
  --markers <markers.json> \
  --groupby <leiden> \
  --output  <annotated.h5ad>
```

**`markers.json` format:**
```json
{
  "T cells": { "pos": ["CD3D", "CD3E", "CD3G"], "neg": [] },
  "B cells": { "pos": ["MS4A1", "CD79A"],        "neg": ["CD3D"] }
}
```

**Scoring algorithm:** For each cluster C and cell type T: `score(C, T) = mean_expr(pos_genes in C) − mean_expr(neg_genes in C)`. Assign the cell type with the highest score if `score > 0`, otherwise label the cluster `"Unknown"`. Store per-cell annotation in `colData$sctype_annotation`.

- [ ] **Step 1: Write the failing test**

```R
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_sctype.R')"
```

Expected: FAIL — `src/scanpy_workflow/annotation/sctype.R` does not exist.

- [ ] **Step 3: Implement `src/scanpy_workflow/annotation/sctype.R`**

```R
#!/usr/bin/env Rscript
# sctype.R — scType-style automated cell type annotation for Nextflow DSL2
#
# Usage:
#   Rscript sctype.R \
#     --input   <clustered.h5ad> \
#     --markers <markers.json> \
#     --groupby <leiden> \
#     --output  <annotated.h5ad>
#
# markers.json format:
#   { "T cells": { "pos": ["CD3D", "CD3E"], "neg": [] }, ... }
#
# Algorithm:
#   For each cluster C and each cell type T:
#     score(C, T) = mean_expr(pos_genes in C) - mean_expr(neg_genes in C)
#   Assign highest-scoring cell type (score > 0) or "Unknown" to each cluster.
#   All cells in the same cluster receive the same annotation.
#   Result stored in colData$sctype_annotation.
#
# Requires: r-jsonlite (conda), zellkonverter (bioconda)
# Optional: remotes::install_github('IanevskiAleksandr/sc-type') for full scType

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(jsonlite)
  library(Matrix)
})

parse_args_sctype <- function() {
  option_list <- list(
    make_option("--input",   type = "character", help = "Input .h5ad path"),
    make_option("--markers", type = "character", help = "Markers JSON path"),
    make_option("--groupby", type = "character", default = "leiden",
                help = "obs column with cluster assignments [default: leiden]"),
    make_option("--output",  type = "character", help = "Output .h5ad path")
  )
  parser <- OptionParser(option_list = option_list)
  opts   <- parse_args(parser)
  if (is.null(opts$input) || is.null(opts$markers) || is.null(opts$output)) {
    print_help(parser)
    stop("--input, --markers, and --output are required", call. = FALSE)
  }
  opts
}

#' Score each cluster against each cell type from the markers list.
#' @param expr_mat  genes × cells expression matrix.
#' @param clusters  Named character vector mapping cell → cluster label.
#' @param markers   Named list: cell_type → list(pos = ..., neg = ...).
#' @return data.frame: cluster, cell_type, score.
score_clusters <- function(expr_mat, clusters, markers) {
  genes_in_mat <- rownames(expr_mat)
  cluster_ids  <- unique(clusters)
  records <- vector("list", length(cluster_ids) * length(markers))
  k <- 1L
  for (cl in cluster_ids) {
    cell_idx <- which(clusters == cl)
    cl_means <- rowMeans(expr_mat[, cell_idx, drop = FALSE])
    for (ct in names(markers)) {
      pos_g     <- intersect(markers[[ct]]$pos, genes_in_mat)
      neg_g     <- intersect(markers[[ct]]$neg, genes_in_mat)
      pos_score <- if (length(pos_g) > 0) mean(cl_means[pos_g]) else 0
      neg_score <- if (length(neg_g) > 0) mean(cl_means[neg_g]) else 0
      records[[k]] <- data.frame(cluster = cl, cell_type = ct,
                                 score = pos_score - neg_score,
                                 stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  do.call(rbind, records)
}

run_sctype <- function(input_path, markers_path, groupby, output_path) {
  message("sctype.R: reading ", input_path)
  sce <- readH5AD(input_path, reader = "R")

  if (!groupby %in% names(colData(sce))) {
    stop("colData does not contain column '", groupby, "'. ",
         "Available: ", paste(names(colData(sce)), collapse = ", "))
  }

  expr_mat <- if ("log_norm" %in% assayNames(sce)) {
    as.matrix(assay(sce, "log_norm"))
  } else {
    message("sctype.R: 'log_norm' assay not found; using 'X'")
    as.matrix(assay(sce, "X"))
  }

  message("sctype.R: loading markers from ", markers_path)
  # simplifyVector = TRUE (default): JSON arrays become atomic character vectors.
  # simplifyVector = FALSE would return list objects that break intersect().
  markers <- fromJSON(markers_path)
  for (ct in names(markers)) {
    # Coerce to character; handles NULL (absent key) and numeric/list edge cases
    if (is.null(markers[[ct]]$pos)) markers[[ct]]$pos <- character(0)
    if (is.null(markers[[ct]]$neg)) markers[[ct]]$neg <- character(0)
    markers[[ct]]$pos <- as.character(markers[[ct]]$pos)
    markers[[ct]]$neg <- as.character(markers[[ct]]$neg)
  }

  clusters <- setNames(as.character(colData(sce)[[groupby]]), colnames(sce))
  message("sctype.R: scoring ", length(unique(clusters)), " clusters against ",
          length(markers), " cell types")

  scores_df <- score_clusters(expr_mat, clusters, markers)

  # Per cluster: pick highest-scoring cell type (or "Unknown" if max <= 0)
  best_df <- do.call(rbind, lapply(split(scores_df, scores_df$cluster), function(df) {
    best_idx <- which.max(df$score)
    label    <- if (df$score[best_idx] > 0) df$cell_type[best_idx] else "Unknown"
    data.frame(cluster = df$cluster[1], annotation = label, stringsAsFactors = FALSE)
  }))

  cluster_to_annot      <- setNames(best_df$annotation, best_df$cluster)
  sce$sctype_annotation <- cluster_to_annot[clusters]

  message("sctype.R: annotation summary:")
  print(table(sce$sctype_annotation))

  message("sctype.R: writing ", output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeH5AD(sce, output_path)
  message("sctype.R: done")
}

main <- function() {
  opts <- parse_args_sctype()
  run_sctype(opts$input, opts$markers, opts$groupby, opts$output)
}

if (!interactive()) main()
```

- [ ] **Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('tests/unit/r/test_sctype.R')"
```

Expected: All 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/annotation/sctype.R tests/unit/r/test_sctype.R
git commit -m "feat: add scType annotation R wrapper with testthat tests"
```

---

## Chunk 6: CI integration

### Task 8: Wire R tests into CI

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Run full R test suite locally before writing CI**

```bash
Rscript -e "testthat::test_dir('tests/unit/r/', reporter = 'progress')"
```

Expected: All 8 tests across 5 test files PASS.

- [ ] **Step 2: Write `.github/workflows/ci.yml`**

```yaml
# .github/workflows/ci.yml
name: CI

on:
  push:
    branches: ["**"]
  pull_request:
    branches: [main]

jobs:
  python-unit:
    name: Python unit tests
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: Install package
        run: pip install -e ".[dev]"
      - name: Run Python unit tests
        run: pytest tests/unit/ -v --tb=short

  r-unit:
    name: R unit tests
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: "4.3"
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          packages: |
            any::optparse
            any::testthat
            any::jsonlite
            any::remotes
            any::SoupX
            bioc::zellkonverter
            bioc::SingleCellExperiment
            bioc::scran
            bioc::scater
            bioc::celda
      - name: Install GitHub-only packages (ALRA, scType)
        run: make install-r-github
      - name: Run R unit tests
        run: Rscript -e "testthat::test_dir('tests/unit/r/', reporter = 'progress')"

  integration:
    name: Integration tests
    runs-on: ubuntu-latest
    if: github.event_name == 'pull_request'
    needs: [python-unit, r-unit]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: Install package
        run: pip install -e ".[dev]"
      - name: Run integration tests
        run: pytest tests/integration/ -v --tb=short
```

- [ ] **Step 3: Verify CI config is valid YAML**

```bash
python -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))" && echo "YAML valid"
```

Expected: `YAML valid`

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add GitHub Actions workflow for Python and R unit tests"
```
