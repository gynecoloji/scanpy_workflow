# scanpy_workflow

An end-to-end single-cell RNA-seq analysis pipeline built on [Scanpy](https://scanpy.readthedocs.io), [Nextflow DSL2](https://www.nextflow.io), and Bioconductor. It takes raw Cell Ranger output and produces clustered, annotated, and differentially expressed results as `.h5ad` files and an HTML report.

![Pipeline diagram](docs/pipeline.svg)

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Pipeline Steps](#pipeline-steps)
- [QC Report](#qc-report)
- [Pathway Scoring](#pathway-scoring)
- [Differential Expression Methods](#differential-expression-methods)
- [Output Structure](#output-structure)
- [Running Tests](#running-tests)
- [Conda Environments](#conda-environments)
- [CI/CD](#cicd)

---

## Features

| Category | Options |
|---|---|
| QC reporting | Per-sample HTML report: violin, histogram, scatter, knee plot, filter waterfall |
| Ambient RNA removal | SoupX, decontX |
| Doublet detection | Scrublet |
| Normalization | Library-size (scanpy), scran |
| Batch correction | Harmony, BBKNN, scVI |
| Imputation | MAGIC, scVI, ALRA |
| Clustering | Leiden, Louvain |
| Cell-type annotation | CellTypist, scType |
| Differential expression | Wilcoxon, t-test, logistic regression, MAST, DESeq2 pseudo-bulk |
| Pathway scoring | AUCell, ULM (decoupler); scanpy score_genes; MSigDB Hallmark/KEGG/Reactome/GO-BP, PROGENy, custom GMT/CSV |

All inter-process communication is via `.h5ad` files. The pipeline can run with conda, Docker, or Singularity.

---

## Requirements

| Tool | Version |
|---|---|
| Nextflow | ≥ 23.10 |
| conda / mamba | any recent |
| Python | 3.11 (managed by conda env) |
| R | 4.3 (managed by conda env) |

Nextflow and conda are the only tools you need on the host machine.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-org>/scanpy_workflow.git
cd scanpy_workflow
```

### 2. Create conda environments

Three environments cover all pipeline steps:

```bash
# Python — scanpy, harmony, BBKNN, scrublet, celltypist, …
conda env create -f envs/env_scanpy.yaml

# Python + PyTorch — scVI-tools
conda env create -f envs/env_scvi.yaml

# R / Bioconductor — scran, SoupX, decontX, MAST, DESeq2, ALRA, scType
conda env create -f envs/env_r.yaml
```

Install GitHub-only R packages (ALRA, scType) after the conda environment is ready:

```bash
conda run -n env_r Rscript -e "remotes::install_github('nalab-stanford/ALRA')"
conda run -n env_r Rscript -e "remotes::install_github('IanevskiAleksandr/sc-type')"
```

### 3. Install the Python package

```bash
conda activate env_scanpy
pip install -e ".[dev]"
```

---

## Quick Start

### Minimal run (two samples, default settings)

```bash
nextflow run workflow/main.nf \
  -params-file params.yaml \
  -profile conda
```

`params.yaml` points to `data/cellranger/` by default. Place your Cell Ranger output directories there:

```
data/cellranger/
├── sample_A/          # filtered_feature_bc_matrix/ (and raw_feature_bc_matrix/ if ambient removal is on)
└── sample_B/
```

Then set the sample names in `params.yaml`:

```yaml
samples:
  - sample_A
  - sample_B
```

### Resume an interrupted run

```bash
nextflow run workflow/main.nf -params-file params.yaml -profile conda -resume
```

### Docker

```bash
nextflow run workflow/main.nf -params-file params.yaml -profile docker
```

### Singularity / Apptainer — pull from Docker Hub at runtime

```bash
nextflow run workflow/main.nf -params-file params.yaml -profile singularity
nextflow run workflow/main.nf -params-file params.yaml -profile apptainer
```

### Singularity / Apptainer — load pre-built `.sif` files (recommended on HPC)

Build the images once on a machine with internet access:

```bash
apptainer pull /shared/containers/scanpy.sif docker://scanpy_workflow/scanpy:latest
apptainer pull /shared/containers/scvi.sif   docker://scanpy_workflow/scvi:latest
apptainer pull /shared/containers/r_env.sif  docker://scanpy_workflow/r:latest
```

Then point the pipeline at that directory with `sif_dir`:

```bash
nextflow run workflow/main.nf \
  -params-file params.yaml \
  -profile apptainer \
  --sif_dir /shared/containers
```

Or set it permanently in `params.yaml`:

```yaml
sif_dir: "/shared/containers"
```

Expected filenames inside `sif_dir`:

| File | Used by |
|---|---|
| `scanpy.sif` | All Python steps |
| `scvi.sif` | scVI batch correction / imputation |
| `r_env.sif` | scran, SoupX, decontX, MAST, DESeq2, ALRA, scType |

---

## Configuration

All options live in `params.yaml`. Pass it with `-params-file params.yaml`. Every key has a default in `workflow/nextflow.config`; `params.yaml` overrides them.

### Full annotated example

```yaml
# --- I/O ---
input_dir: "data/cellranger/"       # parent directory; each sample is a subdirectory
output_dir: "results/"
samples:
  - sample_A
  - sample_B

# --- Ambient RNA removal (Step 2) ---
ambient:
  enabled: true
  method: soupx                      # soupx | decontx

# --- QC filtering thresholds (Step 4) ---
filter:
  min_genes: 200
  max_genes: 6000
  min_counts: 500
  max_counts: 30000
  max_pct_mito: 20                   # percent mitochondrial reads
  max_pct_ribo: 50                   # percent ribosomal reads
  min_cells: 3                       # minimum cells a gene must appear in

# --- Doublet detection (Step 5) ---
doublet:
  enabled: true
  method: scrublet

# --- Normalization (Step 6) ---
normalization:
  method: library_size               # library_size | scran
  target_sum: 10000                  # used only by library_size

# --- Highly variable genes (Steps 7, 9) ---
hvg:
  n_top_genes: 3000                  # per-sample HVG selection
  post_merge_n_top_genes: 2000       # post-merge HVG selection
  flavor: seurat_v3                  # seurat_v3 | seurat | cell_ranger

# --- PCA (Step 11) ---
pca:
  n_comps: 50

# --- Batch correction (Step 12) ---
batch_correction:
  enabled: true
  method: harmony                    # harmony | scvi | bbknn
  batch_key: sample                  # obs column used as batch label

# --- Imputation (Step 13, optional) ---
imputation:
  enabled: false
  method: magic                      # magic | scvi | alra
  evaluate: true                     # compute imputation quality metrics

# --- Clustering (Steps 14-15) ---
clustering:
  algorithm: leiden                  # leiden | louvain
  resolution: 0.5
  n_neighbors: 15
  evaluate: true                     # compute silhouette / Davies-Bouldin scores

# --- Cell-type annotation (Step 16) ---
annotation:
  unsupervised: true                 # rank marker genes per cluster
  automated: true                    # run automated classifier
  method: celltypist                 # celltypist | sctype
  celltypist_model: "Immune_All_Low.pkl"
  sctype_markers: null               # path to markers JSON; required if method: sctype

# --- Differential expression (Step 17) ---
de:
  method: wilcoxon                   # wilcoxon | t-test | logreg | mast | pseudobulk
  groupby: leiden                    # obs column to group cells by
  sample_key: sample                 # obs column with sample IDs (pseudobulk only)
  min_cells: 10                      # min cells per (cluster, sample) pseudo-bulk

# --- Report (Step 18) ---
report: true
```

---

## Pipeline Steps

The pipeline is split into two sub-workflows and a main orchestrator:

### PER_SAMPLE sub-workflow (Steps 1–7, parallelised per sample)

| Step | Process | Method |
|---|---|---|
| 1 | LOAD | Read Cell Ranger MEX or HDF5; store raw counts in `layers["counts"]` |
| 2 | AMBIENT *(optional)* | SoupX or decontX ambient RNA removal |
| 3 | QC | Compute `n_genes`, `n_counts`, `pct_counts_mt`, `pct_counts_ribo` |
| 3b | QC\_REPORT | Per-sample HTML QC report (runs in parallel with FILTER) |
| 4 | FILTER | Hard thresholds on QC metrics |
| 5 | DOUBLETS *(optional)* | Scrublet doublet scores; flag or remove doublets |
| 6 | NORMALIZE | Library-size (`log1p` to `layers["log_norm"]`) or scran |
| 7 | HVG\_PER\_SAMPLE | Select highly variable genes per sample |

### INTEGRATION sub-workflow (Steps 8–13)

| Step | Process | Method |
|---|---|---|
| 8 | MERGE | Concatenate per-sample `.h5ad`; add `batch_key` column |
| 9 | HVG\_POST\_MERGE | Consensus HVG selection across samples |
| 10 | SCALE *(skipped for scVI path)* | Z-score; clip at `max_value=10` |
| 11 | PCA | `n_comps` principal components |
| 12 | BATCH\_CORRECT *(optional)* | Harmony (corrects `X_pca`), BBKNN (graph-level), or scVI (latent space) |
| 13 | IMPUTE *(optional)* | MAGIC, scVI, or ALRA |

### Downstream (Steps 14–18, main.nf)

| Step | Process | Notes |
|---|---|---|
| 14–15 | CLUSTER | KNN graph → Leiden/Louvain; optional silhouette evaluation |
| 16a | RANK\_GENES | Unsupervised marker gene ranking per cluster |
| 16b | ANNOTATE | CellTypist (pre-trained model) or scType (custom markers) |
| 16c | PATHWAY\_SCORE *(optional)* | Per-cell pathway scoring; runs in parallel with 16a/16b |
| 17 | DE | Differential expression (see below) |
| 18 | REPORT | HTML report + final `.h5ad` |

---

## QC Report

A per-sample HTML QC report is generated automatically at step 3b for every sample, using the pre-filter AnnData so the full cell distribution is visible. The report is self-contained (plots embedded as base64 PNG) and requires no web server.

### Contents

| Section | Description |
|---|---|
| Summary cards | Total cells, cells passing filters (%), median genes/cell, median UMI/cell, median % mito |
| Filter waterfall | Bar chart of cells remaining after each sequential filter + table with cells removed per step |
| Violin plots | `n_genes`, UMI counts, % mito, % ribo — with dashed threshold lines |
| Histograms | Same 4 metrics as bar histograms with threshold lines |
| Scatter plot | Total counts vs genes detected, coloured by % mito, threshold lines overlaid |
| Knee plot | Cell rank vs total UMI (log–log) — helps identify empty droplets |
| Metric statistics | Median, mean, min, max, 5th/95th percentiles per metric |
| Thresholds table | Exact filter values applied |

### Output location

```
results/per_sample/<sample>/qc_report/
└── qc_report.html      # self-contained, single-file report
```

### Running standalone

```bash
scanpy-workflow qc-report \
  --input        data/cellranger/sample_A/filtered_feature_bc_matrix/ \
  --sample       sample_A \
  --output       qc_report_A.html \
  --min-genes    200 \
  --max-genes    6000 \
  --min-counts   500 \
  --max-counts   30000 \
  --max-pct-mito 20 \
  --max-pct-ribo 50
```

The `--input` can be a Cell Ranger directory **or** an `.h5ad` file that has already had `compute_qc_metrics()` applied.

---

## Pathway Scoring

Pathway scoring runs per-cell on the clustered `.h5ad` (step 16c) and is independent of annotation and DE. It can be enabled for any combination of built-in or custom gene sets.

### Enabling pathway scoring

```yaml
pathway:
  enabled: true
  source: msigdb_hallmark   # see sources below
  method: aucell            # see methods below
  groupby: leiden
  organism: human           # human | mouse
  min_n: 5                  # drop gene sets with fewer than this many genes in the data
  custom_genesets: null     # path to .gmt or .csv (required when source: custom)
```

### Gene set sources

| `source` | Description | Requires |
|---|---|---|
| `msigdb_hallmark` | MSigDB Hallmark collection (50 gene sets) | decoupler-py |
| `msigdb_kegg` | MSigDB KEGG pathway collection | decoupler-py |
| `msigdb_reactome` | MSigDB Reactome pathway collection | decoupler-py |
| `msigdb_gobp` | MSigDB GO Biological Process collection | decoupler-py |
| `progeny` | PROGENy — 14 cancer signalling pathways with curated weights | decoupler-py |
| `custom` | User-supplied file (`.gmt` or `.csv`/`.tsv`) | — |

### Scoring methods

| `method` | Algorithm | Best for | Requires |
|---|---|---|---|
| `aucell` | Rank-based area under recovery curve | Sparse scRNA-seq; robust to dropouts | decoupler-py |
| `ulm` | Univariate linear model; uses gene weights | Weighted gene sets (PROGENy, DoRothEA) | decoupler-py |
| `scanpy_score` | Mean expression of gene set minus control genes | Quick exploration; no extra dependency | — |

### Custom gene sets

**GMT format** (standard MSigDB export):
```
HALLMARK_APOPTOSIS	Brief description	CASP3	CASP7	CASP9	...
MY_PATHWAY	My custom set	GeneA	GeneB	GeneC
```

**CSV format** (columns `gene_set` and `gene`; optional `weight`):
```csv
gene_set,gene,weight
MY_PATHWAY,GeneA,1.0
MY_PATHWAY,GeneB,0.8
OTHER_PATH,GeneC,1.0
```

```yaml
pathway:
  enabled: true
  source: custom
  method: aucell
  custom_genesets: "data/my_genesets.gmt"
```

### Output

```
results/pathway/
├── pathway_scored.h5ad         # h5ad with obsm["pathway_scores"]
└── pathway_results/
    ├── pathway_scores.csv      # cells × pathways matrix
    └── cluster_pathway_scores.csv  # per-cluster mean scores
```

`adata.obsm["pathway_scores"]` is a `pandas.DataFrame` (cells × pathways) accessible in Python:

```python
import anndata as ad
adata = ad.read_h5ad("results/pathway/pathway_scored.h5ad")
scores = adata.obsm["pathway_scores"]   # DataFrame: cells × pathways
scores.groupby(adata.obs["leiden"]).mean()  # per-cluster activity
```

### Installing the decoupler extra

```bash
pip install "scanpy-workflow[pathway]"
# or via conda (already included in env_scanpy.yaml):
conda env create -f envs/env_scanpy.yaml
```

---

## Differential Expression Methods

All DE methods write a CSV with the same columns:

| Column | Description |
|---|---|
| `group` | Cluster label (one-vs-rest) |
| `gene` | Gene name |
| `score` | Method-specific ranking score |
| `logfoldchange` | Log fold change (in-cluster mean minus rest mean, or log2FC for DESeq2) |
| `pval` | Raw p-value |
| `pval_adj` | BH-adjusted p-value (per cluster) |

### Wilcoxon / t-test / logistic regression

Scanpy's `rank_genes_groups` implementation. Fast, no extra dependencies. Good for large datasets.

```yaml
de:
  method: wilcoxon   # or t-test | logreg
  groupby: leiden
```

### MAST (hurdle model)

Fits a two-component hurdle model (continuous + discrete) per cluster using [MAST](https://bioconductor.org/packages/MAST/). Accounts for zero-inflation and cellular detection rate (`ngeneson` covariate). Score = `sign(LFC) × √λ_hurdle`.

```yaml
de:
  method: mast
  groupby: leiden
```

Requires: `bioconductor-mast` (included in `env_r.yaml`).

### Pseudo-bulk DESeq2

Aggregates raw counts per (cluster × sample), then runs [DESeq2](https://bioconductor.org/packages/DESeq2/) Wald test. Recommended for multi-sample experiments. Clusters with fewer than 2 samples meeting the `min_cells` threshold are skipped. Score = Wald statistic.

```yaml
de:
  method: pseudobulk
  groupby: leiden
  sample_key: sample   # obs column containing sample IDs
  min_cells: 10        # minimum cells per (cluster, sample) pseudo-replicate
```

Requires: `bioconductor-deseq2` (included in `env_r.yaml`) and at least 2 samples in your dataset.

---

## Output Structure

```
results/
├── de/
│   └── de_results.csv               # DE table (all clusters)
├── clustering/
│   ├── clustered.h5ad
│   └── cluster_evaluation/
├── annotation/
│   ├── annotated.h5ad
│   └── celltypist_results/
├── report/
│   └── report.html                  # Summary report
├── nextflow_report.html             # Nextflow execution report
└── timeline.html                    # Nextflow timeline
```

The final `.h5ad` contains:

| Slot | Content |
|---|---|
| `adata.X` | `log1p`-normalised expression (`layers["log_norm"]`) |
| `layers["counts"]` | Raw integer counts |
| `layers["norm"]` | Normalised counts (pre-log) |
| `layers["log_norm"]` | log1p-normalised counts |
| `obsm["X_pca"]` | PCA embedding |
| `obsm["X_umap"]` | UMAP embedding |
| `obs["leiden"]` | Cluster labels |
| `obs["celltypist_cell_type"]` | CellTypist predictions (if enabled) |
| `obs["sctype_cell_type"]` | scType predictions (if enabled) |

---

## Running Tests

### Python unit tests

```bash
conda activate env_scanpy
pip install -e ".[dev]"
pytest tests/unit/ -v
```

Skipped automatically if optional environments (scVI, MAGIC) are not available.

### R unit tests

```bash
conda activate env_r
Rscript -e "testthat::test_dir('tests/unit/r/', reporter = 'progress')"
```

Individual test files:

```bash
Rscript tests/unit/r/test_mast.R
Rscript tests/unit/r/test_pseudobulk.R
Rscript tests/unit/r/test_scran.R
Rscript tests/unit/r/test_soupx.R
Rscript tests/unit/r/test_decontx.R
Rscript tests/unit/r/test_alra.R
Rscript tests/unit/r/test_sctype.R
```

### Integration / smoke test

Requires Nextflow and the `env_scanpy` conda environment:

```bash
conda activate env_scanpy
pytest tests/integration/ -v
```

---

## Conda Environments

| File | Label in Nextflow | Used by |
|---|---|---|
| `envs/env_scanpy.yaml` | `scanpy` | All Python steps (load, QC, normalize, HVG, PCA, batch-correct, cluster, annotate, DE, report) |
| `envs/env_scvi.yaml` | `scvi` | scVI batch correction, scVI imputation |
| `envs/env_r.yaml` | `r_env` | scran, SoupX, decontX, ALRA, scType, MAST, DESeq2 |

Nextflow selects the correct environment automatically based on the `label` declared in each module. You only need the environments for the steps you enable.

---

## Profiles

| Profile | Description |
|---|---|
| `conda` | Uses conda YAML files in `envs/`; recommended for development |
| `docker` | Pulls Docker images; requires Docker daemon |
| `singularity` | Uses Singularity; pulls from Docker Hub unless `sif_dir` is set |
| `apptainer` | Uses Apptainer (successor to Singularity); same image sources and `sif_dir` logic |
| `test` | Minimal CI profile; expects conda envs already activated |

> **Apptainer vs Singularity:** Apptainer is the Linux Foundation fork of Singularity (v1.0+). Both use `.sif` files and the same container images. Use `-profile apptainer` on systems where the binary is `apptainer`, and `-profile singularity` where it is `singularity`. Requires Nextflow ≥ 22.10 for the `apptainer` profile.

---

## CI/CD

GitHub Actions runs three jobs on every push:

| Job | Trigger | What it does |
|---|---|---|
| `python-unit` | every push | `pytest tests/unit/` with Python 3.11 |
| `r-unit` | every push | `testthat::test_dir('tests/unit/r/')` with R 4.3 |
| `integration` | pull request to `main` | Full Nextflow smoke test with conda |

See `.github/workflows/ci.yml` for details.
