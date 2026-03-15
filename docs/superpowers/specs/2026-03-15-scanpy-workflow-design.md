# Scanpy End-to-End Single-Cell Workflow — Design Spec

**Date:** 2026-03-15
**Status:** Approved

---

## 1. Overview

An end-to-end Nextflow pipeline for single-cell RNA-seq analysis starting from 10x Genomics Cell Ranger output. The pipeline supports multiple samples, pluggable batch correction and imputation methods, and produces `.h5ad` AnnData objects, diagnostic plots, and an HTML summary report.

---

## 2. Architecture

**Python package + Nextflow orchestration (Option B).**

- A `scanpy_workflow` Python package contains all analytical logic, organized by module. Each module is independently importable and testable.
- Nextflow orchestrates execution: parallelizes per-sample steps, chains dependent steps, and manages I/O between processes.
- Each Nextflow process calls a CLI entrypoint from the package (e.g., `scanpy-workflow qc --input ... --output ...`).
- A single `params.yaml` config file drives all runtime choices.
- No in-memory state is passed between Nextflow processes — each step reads/writes `.h5ad` files, making the pipeline fully resumable with `-resume`.

### Environment Isolation

Each Nextflow process declares both a `conda` directive and a `container` directive. A single profile flag switches between them:

| Profile | Usage |
|---------|-------|
| `-profile conda` | Local development |
| `-profile docker` | Docker containers |
| `-profile singularity` | HPC/cluster |

| Conda env | Language | Used by |
|-----------|----------|---------|
| `env_scanpy` | Python | QC, filtering, doublets, preprocessing, clustering, annotation, DE, reporting |
| `env_scvi` | Python (torch) | scVI/scANVI batch correction and imputation |
| `env_r` | R | SoupX, DecontX, ALRA, scType |

---

## 3. Pipeline Steps

### Per-Sample (parallelized)

| Step | Tool | Optional |
|------|------|----------|
| 1. Load Cell Ranger output | scanpy `read_10x_h5` / `read_10x_mtx` | No |
| 2. Ambient RNA removal | SoupX (R) or DecontX (R), pluggable | Yes |
| 3. QC metrics | scanpy `pp.calculate_qc_metrics` | No |
| 4. Filter | User-defined cutoffs from `params.yaml` | No |
| 5. Doublet detection | Scrublet | Yes (enabled by default) |
| 6. Normalization + log1p | `pp.normalize_total` + `pp.log1p` | No |
| 7. HVG selection | `pp.highly_variable_genes` | No |
| 8. Scale (per-sample) | `pp.scale` | No |

### Shared (after merge)

| Step | Tool | Optional |
|------|------|----------|
| 9. Merge samples | `ad.concat` | No |
| 10. PCA | `pp.pca` | No |
| 11. Batch correction | Harmony / scVI / BBKNN, pluggable | **Yes** |
| 11a. Batch correction evaluation | kBET + LISI scores + UMAP pre/post | If step 11 enabled |
| 12. Imputation | MAGIC / scVI / ALRA (R), pluggable | **Yes** |
| 12a. Imputation evaluation | Dropout correlation + distribution plots + UMAP comparison | If step 12 enabled |
| 13. Neighborhood graph | `pp.neighbors` | No |
| 14. Leiden/Louvain clustering | `tl.leiden` / `tl.louvain` | No |
| 14a. Clustering evaluation | Silhouette, Davies-Bouldin, Calinski-Harabasz, marker dot plots | No |
| 15a. Unsupervised annotation | `tl.rank_genes_groups` (marker genes) | No |
| 15b. Automated annotation | CellTypist or scType (R), parallel with 15a | No |
| 16. Differential expression | Wilcoxon / t-test / logreg | No |
| 17. Report | HTML report + final `.h5ad` | No |

---

## 4. Evaluation Steps

### Batch Correction Evaluation (step 11a)
- kBET score (batch mixing)
- iLISI / cLISI (integration and cell-type LISI)
- UMAP plots: pre-correction vs post-correction, colored by batch and cell type

### Imputation Evaluation (step 12a)
- Pearson/Spearman correlation on held-out dropout values (mask known zeros, test recovery)
- Gene expression distribution plots (pre vs post imputation)
- UMAP comparison (pre vs post imputation)

### Clustering Evaluation (step 14a)
- Silhouette score
- Davies-Bouldin index
- Calinski-Harabasz index
- Marker gene dot plots per cluster

---

## 5. Configuration

All runtime parameters live in a single `params.yaml`:

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
  method: scrublet

# Normalization
normalization:
  method: library_size           # library_size | scran
  target_sum: 10000

# HVG
hvg:
  n_top_genes: 3000
  flavor: seurat_v3              # seurat_v3 | cell_ranger | seurat

# PCA
pca:
  n_comps: 50

# Batch correction (optional)
batch_correction:
  enabled: true
  method: harmony                # harmony | scvi | bbknn
  batch_key: sample

# Imputation (optional)
imputation:
  enabled: false
  method: magic                  # magic | scvi | alra

# Clustering
clustering:
  algorithm: leiden              # leiden | louvain
  resolution: 0.5
  n_neighbors: 15

# Annotation
annotation:
  unsupervised: true
  automated: true
  method: celltypist             # celltypist | sctype
  celltypist_model: "Immune_All_Low.pkl"

# Differential expression
de:
  method: wilcoxon               # wilcoxon | t-test | logreg
  groupby: leiden

# Output
output_dir: "results/"
report: true
```

---

## 6. Output Structure

```
results/
├── per_sample/{sample}/
│   ├── ambient/                 # decontaminated .h5ad + plots
│   ├── qc/                      # QC metrics + violin/scatter plots
│   ├── filtered/                # filtered .h5ad
│   ├── doublets/                # doublet scores + filtered .h5ad
│   └── preprocessed/           # normalized + HVG + scaled .h5ad
├── merged/                      # merged .h5ad
├── pca/                         # PCA .h5ad + elbow plot
├── batch_correction/            # corrected .h5ad + evaluation plots
├── imputation/                  # imputed .h5ad + evaluation plots
├── clustering/                  # clustered .h5ad + UMAP + evaluation metrics
├── annotation/                  # annotated .h5ad + marker plots
├── de/                          # DE results (CSV) + volcano plots
└── report/
    ├── report.html              # HTML summary report
    └── final.h5ad               # final annotated AnnData
```

---

## 7. Python Package Structure

```
scanpy_workflow/
├── src/
│   └── scanpy_workflow/
│       ├── __init__.py
│       ├── cli.py                  # CLI entrypoints (one per step)
│       ├── io/
│       │   └── load.py             # Load MEX or HDF5 from Cell Ranger
│       ├── qc/
│       │   ├── metrics.py          # Compute QC metrics
│       │   ├── filter.py           # Apply cutoff filters
│       │   └── plots.py            # QC violin/scatter plots
│       ├── ambient/
│       │   ├── soupx.R             # SoupX wrapper
│       │   └── decontx.R           # DecontX wrapper
│       ├── doublets/
│       │   └── scrublet.py
│       ├── preprocessing/
│       │   ├── normalize.py        # Normalization + log1p
│       │   ├── hvg.py              # HVG selection
│       │   ├── scale.py            # Per-sample scaling
│       │   └── pca.py
│       ├── integration/
│       │   ├── harmony.py
│       │   ├── scvi.py
│       │   ├── bbknn.py
│       │   └── evaluate.py         # kBET, LISI, UMAP comparison
│       ├── imputation/
│       │   ├── magic.py
│       │   ├── scvi.py
│       │   ├── alra.R              # ALRA wrapper
│       │   └── evaluate.py         # Dropout correlation, distribution plots
│       ├── clustering/
│       │   ├── cluster.py          # Leiden/Louvain + neighbor graph
│       │   └── evaluate.py         # Silhouette, DB, CH, marker dot plots
│       ├── annotation/
│       │   ├── markers.py          # rank_genes_groups
│       │   ├── celltypist.py
│       │   └── sctype.R            # scType wrapper
│       ├── de/
│       │   └── de.py               # Differential expression
│       ├── reporting/
│       │   └── report.py           # HTML report generation (Jinja2)
│       └── utils/
│           ├── io.py               # AnnData read/write helpers
│           └── config.py           # params.yaml loader
├── workflow/
│   ├── main.nf
│   ├── nextflow.config             # profiles: conda, docker, singularity
│   ├── modules/
│   │   ├── load.nf
│   │   ├── ambient.nf
│   │   ├── qc.nf
│   │   ├── filter.nf
│   │   ├── doublets.nf
│   │   ├── preprocess.nf
│   │   ├── merge.nf
│   │   ├── pca.nf
│   │   ├── batch_correction.nf
│   │   ├── imputation.nf
│   │   ├── clustering.nf
│   │   ├── annotation.nf
│   │   ├── de.nf
│   │   └── report.nf
│   └── subworkflows/
│       ├── per_sample.nf           # steps 1-8 per sample
│       └── integration.nf          # merge → PCA → batch correction → imputation
├── envs/
│   ├── env_scanpy.yaml
│   ├── env_scvi.yaml
│   └── env_r.yaml
├── containers/
│   ├── Dockerfile.scanpy
│   ├── Dockerfile.scvi
│   └── Dockerfile.r
├── tests/
│   ├── unit/
│   └── integration/
├── docs/
└── pyproject.toml
```

---

## 8. Testing Strategy

### Unit Tests (`tests/unit/`)
- One test file per module using small synthetic AnnData objects
- Test each step in isolation (e.g., filter correctly removes cells outside cutoff bounds)
- Test pluggable method dispatch (e.g., `method=harmony` calls the correct function)
- Test R wrappers produce expected output format

### Integration Tests (`tests/integration/`)
- End-to-end run on a small public dataset (PBMC 3k from 10x Genomics)
- Two-sample run to validate merge + batch correction path
- Run with all optional steps enabled and disabled
- Validate final `.h5ad` contains expected `.obs`, `.obsm`, `.uns` keys

### Evaluation Tests
- Assert evaluation metrics are computed and written for each optional evaluation step
- Do not assert correctness of scores — assert presence and format of outputs

### CI
- GitHub Actions: unit tests on every push, integration tests on PRs to `main`
