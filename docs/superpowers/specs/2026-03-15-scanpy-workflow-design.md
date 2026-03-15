# Scanpy End-to-End Single-Cell Workflow — Design Spec

**Date:** 2026-03-15
**Status:** Approved

---

## 1. Overview

An end-to-end Nextflow **DSL2** pipeline for single-cell RNA-seq analysis starting from 10x Genomics Cell Ranger output. The pipeline supports multiple samples, pluggable batch correction and imputation methods, and produces `.h5ad` AnnData objects, diagnostic plots, and an HTML summary report.

---

## 2. Architecture

**Python package + Nextflow DSL2 orchestration (Option B).**

- A `scanpy_workflow` Python package contains all analytical logic, organized by module. Each module is independently importable and testable.
- Nextflow (DSL2) orchestrates execution: parallelizes per-sample steps, chains dependent steps, and manages I/O between processes.
- A single `params.yaml` config file drives all runtime choices.
- No in-memory state is passed between Nextflow processes — each step reads/writes `.h5ad` files, making the pipeline fully resumable with `-resume`.
- On startup, `config.py` validates `params.yaml` for invalid combinations (e.g., scVI used for both batch correction and imputation). Checks that reference `.obs` columns created by later steps (e.g., `de.groupby`) are **warning-only** — they cannot be validated at startup since the column does not yet exist.

### Nextflow–Python/R Dispatch Pattern

- **Python steps:** Each Nextflow process calls a Python CLI entrypoint (e.g., `scanpy-workflow qc --input ... --output ...`) registered in `cli.py`.
- **R steps:** Nextflow processes for R-backed steps call `Rscript path/to/script.R --input ... --output ...` **directly** — they do NOT route through `cli.py`. The `ambient`, `normalize` (scran), `impute` (ALRA), and `annotate` (scType) modules each declare `env_r` and invoke the R script directly.
- **Dual-dispatch modules:** When a step supports both Python and R backends (e.g., `annotation.nf` supports CellTypist via Python CLI and scType via `Rscript`), the `.nf` module contains an `if/else` branch on `params.annotation.method`:
  ```
  if (params.annotation.method == "celltypist") {
      // call scanpy-workflow annotate
  } else if (params.annotation.method == "sctype") {
      // call Rscript sctype.R
  }
  ```
  The same pattern applies to `imputation.nf` (MAGIC/scVI via Python; ALRA via Rscript) and `ambient.nf` (SoupX vs DecontX, both R).

### Environment Isolation

Each Nextflow process declares both a `conda` directive and a `container` directive. A single profile flag switches between them:

| Profile | Usage |
|---------|-------|
| `-profile conda` | Local development |
| `-profile docker` | Docker containers |
| `-profile singularity` | HPC/cluster |

| Conda env | Language | Used by | Key packages |
|-----------|----------|---------|--------------|
| `env_scanpy` | Python | QC, filtering, doublets, normalization (library_size), HVG, scale, PCA, Harmony, BBKNN, MAGIC, clustering, CellTypist, rank_genes_groups, DE, reporting, batch correction evaluation | `scanpy`, `anndata`, `harmony-pytorch`, `bbknn`, `magic-impute`, `scrublet`, `celltypist`, `scib-metrics` |
| `env_scvi` | Python (torch) | scVI/scANVI batch correction and imputation | `scvi-tools`, `torch` |
| `env_r` | R | SoupX, DecontX (via celda), scran, ALRA, scType | `SoupX`, `bioconductor-celda`, `scran`, `ALRA`, `scType`, `zellkonverter` |

---

## 3. Pipeline Steps

### Per-Sample (parallelized)

| Step | Tool | Env | Optional |
|------|------|-----|----------|
| 1. Load Cell Ranger output | scanpy `read_10x_h5` / `read_10x_mtx` | `env_scanpy` | No |
| 2. Ambient RNA removal | SoupX (R) or DecontX via `bioconductor-celda` (R), pluggable | `env_r` | Yes |
| 3. QC metrics | scanpy `pp.calculate_qc_metrics` | `env_scanpy` | No |
| 4. Filter | User-defined cutoffs from `params.yaml` | `env_scanpy` | No |
| 5. Doublet detection | Scrublet (sole supported method; `doublet.method` reserved for future expansion) | `env_scanpy` | Yes (default: true) |
| 6. Normalization + log1p | `pp.normalize_total` + `pp.log1p` (library_size) or `scran.R` (scran) | `env_scanpy` / `env_r` | No |
| 7. HVG selection (per-sample) | `pp.highly_variable_genes` | `env_scanpy` | No |

> **Note on layer preservation (library_size path):** After step 6:
> - `layers["counts"]` — raw integer counts
> - `layers["norm"]` — library-size normalized (post `normalize_total`, pre log1p)
> - `layers["log_norm"]` — log1p-transformed normalized counts
> - `adata.X` set to `layers["log_norm"]`

> **Note on layer preservation (scran path):** `scran.R` must write the same three layers to maintain a consistent contract with downstream Python steps: `layers["counts"]` (raw), `layers["norm"]` (scran size-factor normalized, before log), `layers["log_norm"]` (log1p of scran-normalized), and set `adata.X = layers["log_norm"]`. This contract ensures `pca.py`, `scale.py`, and all downstream Python steps can assume the same layer structure regardless of normalization method.

> **Note on HVG flavor and layer selection:** When `hvg.flavor: seurat_v3`, `hvg.py` passes `layer="counts"` to `pp.highly_variable_genes`. When `hvg.flavor: seurat` or `cell_ranger`, `hvg.py` uses `adata.X` (log-normalized).

### Shared (after merge)

| Step | Tool | Env | Optional |
|------|------|-----|----------|
| 8. Merge samples | `ad.concat` | `env_scanpy` | No |
| 9. HVG selection (post-merge) | `pp.highly_variable_genes` | `env_scanpy` | No |
| 10. Scale | `pp.scale`; **skipped if scVI batch correction** | `env_scanpy` | No |
| 11. PCA | `pp.pca` (see note on `adata.X` contract) | `env_scanpy` | No |
| 12. Batch correction + evaluation | Harmony / scVI / BBKNN, pluggable; evaluation embedded in same process | `env_scanpy` / `env_scvi` | **Yes** |
| 13. Imputation + evaluation | MAGIC / scVI / ALRA (R), pluggable; evaluation embedded in same process (optional) | `env_scanpy` / `env_scvi` / `env_r` | **Yes** |
| 14. Neighborhood graph | `pp.neighbors` using `use_rep` from `adata.uns["neighbors_use_rep"]`; skipped for BBKNN | `env_scanpy` | No |
| 15. Leiden/Louvain clustering | `tl.leiden` / `tl.louvain` | `env_scanpy` | No |
| 15a. Clustering evaluation | Silhouette, Davies-Bouldin, Calinski-Harabasz, marker dot plots | `env_scanpy` | Yes (default: true) |
| 16a. Unsupervised annotation | `tl.rank_genes_groups` (marker genes) | `env_scanpy` | No |
| 16b. Automated annotation | CellTypist (Python) or scType (R), parallel with 16a | `env_scanpy` / `env_r` | No |
| 17. Differential expression | Wilcoxon / t-test / logreg | `env_scanpy` | No |
| 18. Report | HTML report + final `.h5ad` | `env_scanpy` | No |

> **Note on Merge and `batch_key`:** `merge.py` calls `ad.concat(adatas, label=params.batch_correction.batch_key, keys=sample_names)`, creating `adata.obs[batch_key]` (default: `"sample"`). The `batch_key` in `params.yaml` defines both the obs column name written here and the key used by the batch correction method.

> **Note on PCA `adata.X` contract:**
> - Non-scVI path: `scale.py` sets `adata.X` to the scaled data. `pca.py` calls `pp.pca(adata, use_highly_variable=True)` on this scaled `adata.X`.
> - scVI path (scaling skipped): `adata.X` remains `layers["log_norm"]`. `pca.py` explicitly verifies (assert) that `adata.X` equals `adata.layers["log_norm"]` before calling `pp.pca`.

> **Note on `use_rep` handoff:** Each batch correction step writes `adata.uns["neighbors_use_rep"]` before saving. `cluster.py` reads this key to determine `use_rep` for `pp.neighbors`.
>
> | Method | Input embedding | `adata.uns["neighbors_use_rep"]` | `pp.neighbors` behavior |
> |--------|----------------|----------------------------------|------------------------|
> | Harmony | `X_pca` | `"X_pca_harmony"` | Called with `use_rep="X_pca_harmony"` |
> | scVI | `layers["counts"]` (raw counts required by scVI VAE) | `"X_scVI"` | Called with `use_rep="X_scVI"` |
> | BBKNN | `X_pca` (called as `bbknn.bbknn(adata, use_rep="X_pca")`) | `"skip"` | Skipped — BBKNN writes `obsp["connectivities"]` + `obsp["distances"]` directly |
> | None | — | `"X_pca"` | Called with `use_rep="X_pca"` |

> **Note on scVI for batch correction and imputation:**
> - `batch_correction.method: scvi` and `imputation.method: scvi` **cannot both be set** — `config.py` raises a hard error at startup.
> - When `batch_correction.method: scvi`, scVI denoised expression is stored in `adata.layers["scvi_denoised"]` automatically.
> - When `imputation.method: scvi` (and `batch_correction.method` is NOT `scvi`), `imputation/scvi.py` trains its own model using `layers["counts"]` (raw counts required by scVI VAE), stores denoised output in `layers["scvi_denoised"]`, and saves the model to `results/imputation/scvi_model/`.

> **Note on MAGIC input layer:** MAGIC operates on `adata.layers["norm"]` (library-size normalized, pre-log) per MAGIC's recommended usage.

> **Note on `de.groupby`:** If not set, defaults to `clustering.algorithm` at runtime. `config.py` emits a **warning** (not a hard error) if `de.groupby` does not exist as an `.obs` column at startup — it will be created by the clustering step.

> **Note on clustering evaluation performance:** Silhouette score and Davies-Bouldin are O(n²). `evaluate.py` subsamples to ≤20,000 cells (`random_state=42`) for these two metrics. Calinski-Harabasz runs on the full dataset.

---

## 4. Evaluation Steps

### Batch Correction Evaluation (step 12, embedded, always runs when BC is enabled)
Evaluation is embedded in the same Nextflow process as batch correction. No separate toggle.
- Uses `scib-metrics` (Python, `env_scanpy`): kBET approximation, iLISI, cLISI, ASW batch, graph connectivity
- UMAP plots: pre-correction vs post-correction, colored by batch and cell type

### Imputation Evaluation (step 13, embedded, controlled by `imputation.evaluate`)
Evaluation is embedded in the same Nextflow process as imputation, and runs only when `imputation.evaluate: true`.
- **In-silico dropout evaluation:** a deep copy of `adata.layers["norm"]` is made; 10% of non-zero values are randomly masked to zero (`random_state=42`); imputation runs on the copy; Pearson and Spearman correlation are computed between masked true values and imputed values at those positions. The production `.h5ad` uses unmasked data.
- Gene expression distribution plots (pre vs post)
- UMAP comparison (pre vs post)

### Clustering Evaluation (step 15a, controlled by `clustering.evaluate`)
- Silhouette score + Davies-Bouldin: subsampled to ≤20,000 cells, `random_state=42`
- Calinski-Harabasz: full dataset
- Marker gene dot plots per cluster

---

## 5. Configuration

```yaml
# Input
input_dir: "data/cellranger/"
samples:
  - sample_A
  - sample_B

# Ambient RNA removal (optional)
ambient:
  enabled: true
  method: soupx                  # soupx | decontx

# QC & Filter
filter:
  min_genes: 200
  max_genes: 6000
  min_counts: 500
  max_counts: 30000
  max_pct_mito: 20
  max_pct_ribo: 50
  min_cells: 3

# Doublet detection
doublet:
  enabled: true
  method: scrublet               # scrublet (reserved; sole supported method)

# Normalization
normalization:
  method: library_size           # library_size | scran
  target_sum: 10000              # library_size only; scran ignores this

# HVG
hvg:
  n_top_genes: 3000              # per-sample
  post_merge_n_top_genes: 2000   # post-merge
  flavor: seurat_v3              # seurat_v3 (uses layers["counts"]) | cell_ranger | seurat

# PCA
pca:
  n_comps: 50

# Batch correction (optional)
batch_correction:
  enabled: true
  method: harmony                # harmony | scvi | bbknn
  batch_key: sample              # obs column created by merge step

# Imputation (optional)
# Note: method: scvi is invalid if batch_correction.method is also scvi
imputation:
  enabled: false
  method: magic                  # magic | scvi | alra
  evaluate: true                 # run in-silico dropout evaluation when imputation is enabled

# Clustering
clustering:
  algorithm: leiden              # leiden | louvain
  resolution: 0.5
  n_neighbors: 15
  evaluate: true                 # silhouette/DB subsampled ≤20k cells

# Annotation
annotation:
  unsupervised: true             # rank_genes_groups (step 16a)
  automated: true                # CellTypist or scType (step 16b)
  method: celltypist             # celltypist | sctype
  celltypist_model: "Immune_All_Low.pkl"

# Differential expression
# de.groupby defaults to clustering.algorithm; startup check is warning-only
de:
  method: wilcoxon               # wilcoxon | t-test | logreg
  groupby: leiden

# Output
output_dir: "results/"
report: true
```

---

## 6. Output Structure

> `scaled/` is omitted when `batch_correction.method: scvi`. `pca/` always exists.

```
results/
├── per_sample/{sample}/
│   ├── ambient/                 # decontaminated .h5ad + plots
│   ├── qc/                      # QC metrics + violin/scatter plots
│   ├── filtered/                # filtered .h5ad
│   ├── doublets/                # doublet scores + filtered .h5ad
│   └── preprocessed/           # normalized .h5ad (layers: counts, norm, log_norm) + HVG flags
├── merged/                      # merged .h5ad; obs[batch_key] column added by ad.concat
├── scaled/                      # scaled .h5ad (omitted if scVI batch correction)
├── pca/                         # PCA .h5ad + elbow plot
├── batch_correction/            # corrected .h5ad + evaluation plots/metrics
├── imputation/                  # imputed .h5ad + evaluation plots (if imputation.evaluate: true)
│   └── scvi_model/              # scVI model artifact (only if imputation.method: scvi)
├── clustering/                  # clustered .h5ad + UMAP + evaluation metrics (if clustering.evaluate: true)
├── annotation/                  # annotated .h5ad + marker dot plots + automated annotation results
├── de/                          # DE CSV + volcano plots (wilcoxon/t-test) or bar plots (logreg)
└── report/
    ├── report.html
    └── final.h5ad
```

---

## 7. Python Package Structure

```
scanpy_workflow/
├── src/
│   └── scanpy_workflow/
│       ├── __init__.py
│       ├── cli.py                  # Python CLI entrypoints; all multi-word commands use hyphens (Click)
│       ├── io/
│       │   └── load.py
│       ├── qc/
│       │   ├── metrics.py
│       │   ├── filter.py
│       │   └── plots.py
│       ├── ambient/
│       │   ├── soupx.R             # Called directly by Nextflow via Rscript
│       │   └── decontx.R           # Uses bioconductor-celda; called directly by Nextflow
│       ├── doublets/
│       │   └── scrublet.py
│       ├── preprocessing/
│       │   ├── normalize.py        # library_size: stores layers["counts"], ["norm"], ["log_norm"];
│       │   │                       # sets adata.X = layers["log_norm"]
│       │   ├── scran.R             # Must write same three layers + set adata.X = layers["log_norm"];
│       │   │                       # called directly by Nextflow via Rscript
│       │   ├── hvg.py              # --mode per-sample|post-merge; seurat_v3: layer="counts"
│       │   ├── scale.py            # Sets adata.X = scaled data
│       │   └── pca.py              # Non-scVI: adata.X is scaled; scVI: asserts adata.X == layers["log_norm"]
│       ├── integration/
│       │   ├── harmony.py          # Input: X_pca; output: obsm["X_pca_harmony"];
│       │   │                       # writes uns["neighbors_use_rep"]="X_pca_harmony";
│       │   │                       # calls evaluate.py inline
│       │   ├── scvi.py             # Input: layers["counts"] (raw counts required by scVI VAE); output: obsm["X_scVI"] + layers["scvi_denoised"];
│       │   │                       # uses VAEModel from utils/scvi_model.py;
│       │   │                       # writes uns["neighbors_use_rep"]="X_scVI";
│       │   │                       # calls evaluate.py inline
│       │   ├── bbknn.py            # Input: X_pca (calls bbknn.bbknn(adata, use_rep="X_pca"));
│       │   │                       # writes obsp["connectivities"] + obsp["distances"] directly;
│       │   │                       # writes uns["neighbors_use_rep"]="skip";
│       │   │                       # calls evaluate.py inline
│       │   └── evaluate.py         # scib-metrics: kBET approx, iLISI, cLISI, ASW, UMAP comparison
│       ├── imputation/
│       │   ├── magic.py            # Input: layers["norm"]; calls evaluate.py inline if imputation.evaluate
│       │   ├── scvi.py             # Input: layers["counts"] (raw counts required by scVI VAE); uses VAEModel from utils/scvi_model.py;
│       │   │                       # stores layers["scvi_denoised"]; saves model artifact;
│       │   │                       # calls evaluate.py inline if imputation.evaluate
│       │   ├── alra.R              # Called directly by Nextflow via Rscript
│       │   └── evaluate.py         # In-silico dropout on deep copy of layers["norm"],
│       │                           # 10% non-zero mask, random_state=42; production .h5ad is unmasked
│       ├── clustering/
│       │   ├── cluster.py          # Reads uns["neighbors_use_rep"]; skips pp.neighbors if "skip";
│       │   │                       # runs leiden/louvain
│       │   └── evaluate.py         # Runs only if clustering.evaluate=true;
│       │                           # Silhouette+DB subsampled ≤20k cells (random_state=42);
│       │                           # CH on full dataset; marker dot plots
│       ├── annotation/
│       │   ├── markers.py          # rank_genes_groups (step 16a); CLI: scanpy-workflow rank-genes
│       │   ├── celltypist.py       # CellTypist (step 16b, Python); CLI: scanpy-workflow annotate
│       │   └── sctype.R            # scType (step 16b, R); called directly by Nextflow via Rscript
│       ├── de/
│       │   └── de.py               # Outputs CSV: gene, score, logfoldchange, pval, pval_adj;
│       │                           # logreg: pval/pval_adj = NaN; report uses bar plots for logreg
│       ├── reporting/
│       │   └── report.py           # HTML report (Jinja2)
│       └── utils/
│           ├── io.py
│           ├── config.py           # params.yaml loader + startup validation (hard errors for invalid
│           │                       # combinations; warning-only for obs columns created later)
│           └── scvi_model.py       # Shared VAEModel class used by integration/scvi.py
│                                   # and imputation/scvi.py
├── workflow/
│   ├── main.nf                     # Nextflow DSL2. Orchestration:
│   │                               #   1. per_sample subworkflow (fan-out, parallel)
│   │                               #   2. integration subworkflow (merge→HVG→scale→PCA→BC→impute)
│   │                               #   3. clustering.nf
│   │                               #   4. annotation.nf (rank_genes + 16b in parallel)
│   │                               #   5. de.nf
│   │                               #   6. report.nf
│   ├── nextflow.config             # Nextflow DSL2; profiles: conda, docker, singularity
│   ├── modules/
│   │   ├── load.nf
│   │   ├── ambient.nf              # if/else on params.ambient.method → Rscript soupx.R or decontx.R
│   │   ├── qc.nf
│   │   ├── filter.nf
│   │   ├── doublets.nf
│   │   ├── preprocess.nf           # normalize (Python or Rscript) + hvg per-sample
│   │   ├── merge.nf
│   │   ├── hvg_postmerge.nf
│   │   ├── scale.nf
│   │   ├── pca.nf
│   │   ├── batch_correction.nf     # if/else on method → harmony | scvi | bbknn
│   │   ├── imputation.nf           # if/else on method → magic (Python) | scvi (Python) | alra (Rscript)
│   │   ├── clustering.nf
│   │   ├── annotation.nf           # Dispatches: rank_genes (Python CLI) + CellTypist (Python CLI) or
│   │   │                           # scType (Rscript) in parallel; if/else on params.annotation.method
│   │   ├── de.nf
│   │   └── report.nf
│   └── subworkflows/
│       ├── per_sample.nf           # steps 1-7 per sample
│       └── integration.nf          # merge → HVG → scale → PCA → batch correction → imputation
├── envs/
│   ├── env_scanpy.yaml             # scanpy, anndata, harmony-pytorch, bbknn, magic-impute,
│   │                               # scrublet, celltypist, scib-metrics
│   ├── env_scvi.yaml               # scvi-tools, torch
│   └── env_r.yaml                  # SoupX, bioconductor-celda, scran, ALRA, scType, zellkonverter
├── containers/
│   ├── Dockerfile.scanpy
│   ├── Dockerfile.scvi
│   └── Dockerfile.r
├── tests/
│   ├── unit/
│   │   ├── (Python test files, one per module)
│   │   └── r/                      # testthat test files for R wrappers
│   └── integration/
│       └── (end-to-end tests using synthetic 2-sample dataset from conftest.py)
├── docs/
└── pyproject.toml
```

### Python CLI Entrypoints

| CLI command | Step | Conda env |
|-------------|------|-----------|
| `scanpy-workflow load` | 1. Load | `env_scanpy` |
| `scanpy-workflow qc` | 3. QC metrics | `env_scanpy` |
| `scanpy-workflow filter` | 4. Filter | `env_scanpy` |
| `scanpy-workflow doublets` | 5. Doublet detection | `env_scanpy` |
| `scanpy-workflow normalize` | 6. Normalization (library_size) | `env_scanpy` |
| `scanpy-workflow hvg --mode per-sample` | 7. Per-sample HVG | `env_scanpy` |
| `scanpy-workflow merge` | 8. Merge; adds obs[batch_key] | `env_scanpy` |
| `scanpy-workflow hvg --mode post-merge` | 9. Post-merge HVG | `env_scanpy` |
| `scanpy-workflow scale` | 10. Scale | `env_scanpy` |
| `scanpy-workflow pca` | 11. PCA | `env_scanpy` |
| `scanpy-workflow batch-correct` | 12. Batch correction + evaluation | `env_scanpy` / `env_scvi` |
| `scanpy-workflow impute` | 13. Imputation + evaluation (Python methods) | `env_scanpy` / `env_scvi` |
| `scanpy-workflow cluster` | 14/15. Neighbors + clustering + evaluation | `env_scanpy` |
| `scanpy-workflow rank-genes` | 16a. Unsupervised annotation | `env_scanpy` |
| `scanpy-workflow annotate` | 16b. Automated annotation (CellTypist) | `env_scanpy` |
| `scanpy-workflow de` | 17. DE | `env_scanpy` |
| `scanpy-workflow report` | 18. Report | `env_scanpy` |

### DE Output Schema

All DE methods write a CSV with columns: `group`, `gene`, `score`, `logfoldchange`, `pval`, `pval_adj`. `group` identifies the cluster each result belongs to (value of `de.groupby`). For `logreg`, `pval` and `pval_adj` are `NaN`. The reporting step generates bar plots of top-scoring genes for `logreg` instead of volcano plots.

---

## 8. Testing Strategy

### Unit Tests (`tests/unit/`)
- One test file per Python module; small synthetic AnnData objects
- Test each step in isolation
- Test config validation catches: scVI conflict, hard-error cases vs. warning-only cases
- Test `uns["neighbors_use_rep"]` written correctly by each batch correction method
- Test `pca.py` `adata.X` assertion for both scVI and non-scVI paths
- Test `de.py` always outputs all required CSV columns (with NaN for logreg)

### R Wrapper Tests (`tests/unit/r/`)
- `testthat` files in `tests/unit/r/`, one per R wrapper
- Each test: reads `.h5ad` via `zellkonverter`, runs wrapper on small synthetic data, validates output
- Run via: `Rscript -e "testthat::test_dir('tests/unit/r/')"`
- Included in CI

### Integration Tests (`tests/integration/`)
- Uses a **programmatically generated synthetic 2-sample dataset** (~500 cells/sample) from `conftest.py` — no external downloads in CI
- Full pipeline run with all options enabled and disabled
- Validate final `.h5ad` keys, output directory structure, and evaluation outputs
- Optional local test against PBMC 3k (downloaded/cached separately via Makefile)

### Evaluation Tests
- Assert evaluation metrics are present and correctly formatted after each evaluation step
- Assert dropout evaluation uses a deep copy (production `.h5ad` is uncontaminated)

### CI
- GitHub Actions (Nextflow DSL2): Python + R unit tests on every push; integration tests on PRs to `main`
