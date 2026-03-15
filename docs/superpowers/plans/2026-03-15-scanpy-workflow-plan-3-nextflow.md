# Scanpy Workflow — Plan 3: Nextflow DSL2 Pipeline

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Nextflow DSL2 pipeline (`workflow/`), conda environment YAMLs (`envs/`), Dockerfiles (`containers/`), and GitHub Actions CI that wires the Python package (Plan 1) and R wrappers (Plan 2) into a complete end-to-end single-cell analysis pipeline.

**Architecture:** Nextflow DSL2 with modules under `workflow/modules/` and two subworkflows (`per_sample.nf`, `integration.nf`). `main.nf` calls subworkflows then chains clustering → annotation → DE → report. Each process declares a `conda` and `container` directive, with profile selection via `-profile conda|docker|singularity`. All inter-process I/O is `.h5ad` files — no in-memory state crosses process boundaries.

**Tech Stack:** Nextflow DSL2 ≥ 23.10, Python package from Plan 1, R wrappers from Plan 2, conda (env_scanpy / env_scvi / env_r), Docker/Singularity containers, params.yaml config

**Spec:** `docs/superpowers/specs/2026-03-15-scanpy-workflow-design.md`

---

## Chunk 1: Config + scaffold

### Task 1: nextflow.config + params.yaml

**Files:**
- Create: `workflow/nextflow.config`
- Create: `params.yaml`

- [ ] **Step 1: Write `workflow/nextflow.config`**

```groovy
// workflow/nextflow.config
// DSL2 is the default in Nextflow >= 22.03; no explicit enable needed.

// Default params — override with params.yaml (-params-file params.yaml)
params {
    input_dir  = "data/cellranger/"
    samples    = []
    output_dir = "results/"

    ambient      = [enabled: false, method: "soupx"]
    filter       = [min_genes: 200, max_genes: 6000, min_counts: 500,
                    max_counts: 30000, max_pct_mito: 20.0, max_pct_ribo: 50.0,
                    min_cells: 3]
    doublet      = [enabled: true,  method: "scrublet"]
    normalization = [method: "library_size", target_sum: 10000]
    hvg          = [n_top_genes: 3000, post_merge_n_top_genes: 2000,
                    flavor: "seurat_v3"]
    pca          = [n_comps: 50]
    batch_correction = [enabled: true,  method: "harmony", batch_key: "sample"]
    imputation   = [enabled: false, method: "magic", evaluate: true]
    clustering   = [algorithm: "leiden", resolution: 0.5, n_neighbors: 15,
                    evaluate: true]
    annotation   = [unsupervised: true, automated: true,
                    method: "celltypist",
                    celltypist_model: "Immune_All_Low.pkl",
                    sctype_markers: null]
    de           = [method: "wilcoxon", groupby: "leiden"]
    report       = true
}

// Per-environment conda YAML paths
def condaDir = "${projectDir}/../envs"

profiles {
    conda {
        conda.enabled = true
        process {
            withLabel: 'scanpy' { conda = "${condaDir}/env_scanpy.yaml" }
            withLabel: 'scvi'   { conda = "${condaDir}/env_scvi.yaml"   }
            withLabel: 'r_env'  { conda = "${condaDir}/env_r.yaml"      }
        }
    }
    docker {
        docker.enabled = true
        process {
            withLabel: 'scanpy' { container = 'scanpy_workflow/scanpy:latest' }
            withLabel: 'scvi'   { container = 'scanpy_workflow/scvi:latest'   }
            withLabel: 'r_env'  { container = 'scanpy_workflow/r:latest'      }
        }
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        process {
            withLabel: 'scanpy' { container = 'docker://scanpy_workflow/scanpy:latest' }
            withLabel: 'scvi'   { container = 'docker://scanpy_workflow/scvi:latest'   }
            withLabel: 'r_env'  { container = 'docker://scanpy_workflow/r:latest'      }
        }
    }
    test {
        // Minimal profile for CI; uses already-activated conda envs
        params.output_dir = "test_results/"
    }
}

// Capture run info
report {
    enabled   = true
    file      = "${params.output_dir}/nextflow_report.html"
}
timeline {
    enabled = true
    file    = "${params.output_dir}/timeline.html"
}
```

- [ ] **Step 2: Write `params.yaml`**

```yaml
# params.yaml — runtime configuration for scanpy_workflow Nextflow pipeline
# Override defaults in nextflow.config. Pass with: nextflow run ... -params-file params.yaml

input_dir: "data/cellranger/"
samples:
  - sample_A
  - sample_B

ambient:
  enabled: true
  method: soupx           # soupx | decontx

filter:
  min_genes: 200
  max_genes: 6000
  min_counts: 500
  max_counts: 30000
  max_pct_mito: 20
  max_pct_ribo: 50
  min_cells: 3

doublet:
  enabled: true
  method: scrublet

normalization:
  method: library_size    # library_size | scran
  target_sum: 10000

hvg:
  n_top_genes: 3000
  post_merge_n_top_genes: 2000
  flavor: seurat_v3       # seurat_v3 | seurat | cell_ranger

pca:
  n_comps: 50

batch_correction:
  enabled: true
  method: harmony         # harmony | scvi | bbknn
  batch_key: sample

imputation:
  enabled: false
  method: magic           # magic | scvi | alra
  evaluate: true

clustering:
  algorithm: leiden       # leiden | louvain
  resolution: 0.5
  n_neighbors: 15
  evaluate: true

annotation:
  unsupervised: true
  automated: true
  method: celltypist      # celltypist | sctype
  celltypist_model: "Immune_All_Low.pkl"
  sctype_markers: null    # path to markers.json; required if method: sctype

de:
  method: wilcoxon        # wilcoxon | t-test | logreg
  groupby: leiden

output_dir: "results/"
report: true
```

- [ ] **Step 3: Commit**

```bash
git add workflow/nextflow.config params.yaml
git commit -m "feat: add nextflow.config and params.yaml scaffold"
```

---

### Task 2: Conda environment YAMLs and Dockerfiles

**Files:**
- Create: `envs/env_scanpy.yaml`
- Create: `envs/env_scvi.yaml`
- Create: `containers/Dockerfile.scanpy`
- Create: `containers/Dockerfile.scvi`
- Create: `containers/Dockerfile.r`

(Note: `envs/env_r.yaml` already created in Plan 2 Task 1.)

- [ ] **Step 1: Write `envs/env_scanpy.yaml`**

```yaml
# envs/env_scanpy.yaml
name: env_scanpy
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - scanpy>=1.9
  - anndata>=0.9
  - harmonypy>=0.0.9
  - bbknn>=1.6
  - magic-impute>=3.0
  - scrublet>=0.2
  - celltypist>=1.5
  - scib-metrics>=0.4
  - leidenalg>=0.10
  - louvain>=0.8
  - python-igraph>=0.11
  - scikit-learn>=1.3
  - click>=8.0
  - pyyaml>=6.0
  - jinja2>=3.0
  - numpy>=1.24
  - pandas>=2.0
  - scipy>=1.10
  - matplotlib>=3.7
  - seaborn>=0.12
  - pip
  - pip:
    - scanpy_workflow  # install local package (editable install handled by CI)
```

- [ ] **Step 2: Write `envs/env_scvi.yaml`**

```yaml
# envs/env_scvi.yaml
name: env_scvi
channels:
  - conda-forge
  - pytorch
  - defaults
dependencies:
  - python=3.11
  - pytorch>=2.0
  - scvi-tools>=1.0
  - anndata>=0.9
  - numpy>=1.24
  - pandas>=2.0
  - click>=8.0
  - pyyaml>=6.0
  - pip
  - pip:
    - scanpy_workflow
```

- [ ] **Step 3: Write `containers/Dockerfile.scanpy`**

```dockerfile
# containers/Dockerfile.scanpy
FROM continuumio/miniconda3:24.1.2-0

WORKDIR /app

COPY envs/env_scanpy.yaml /tmp/env_scanpy.yaml
RUN conda env create -f /tmp/env_scanpy.yaml && conda clean -afy

COPY . /app
RUN conda run -n env_scanpy pip install -e /app --no-deps

SHELL ["conda", "run", "-n", "env_scanpy", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "-n", "env_scanpy", "scanpy-workflow"]
```

- [ ] **Step 4: Write `containers/Dockerfile.scvi`**

```dockerfile
# containers/Dockerfile.scvi
FROM continuumio/miniconda3:24.1.2-0

WORKDIR /app

COPY envs/env_scvi.yaml /tmp/env_scvi.yaml
RUN conda env create -f /tmp/env_scvi.yaml && conda clean -afy

COPY . /app
RUN conda run -n env_scvi pip install -e /app --no-deps

SHELL ["conda", "run", "-n", "env_scvi", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "-n", "env_scvi", "scanpy-workflow"]
```

- [ ] **Step 5: Write `containers/Dockerfile.r`**

```dockerfile
# containers/Dockerfile.r
FROM bioconductor/bioconductor_docker:RELEASE_3_18

WORKDIR /app

COPY envs/env_r.yaml /tmp/env_r.yaml

# Install conda for R env management
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda env create -f /tmp/env_r.yaml && \
    /opt/conda/bin/conda clean -afy

# Install GitHub-only R packages
RUN /opt/conda/envs/env_r/bin/Rscript \
    -e "remotes::install_github('nalab-stanford/ALRA', upgrade='never')" \
    -e "remotes::install_github('IanevskiAleksandr/sc-type', upgrade='never')"

COPY src/scanpy_workflow/ambient/   /app/ambient/
COPY src/scanpy_workflow/preprocessing/ /app/preprocessing/
COPY src/scanpy_workflow/imputation/    /app/imputation/
COPY src/scanpy_workflow/annotation/    /app/annotation/
```

- [ ] **Step 6: Commit**

```bash
git add envs/env_scanpy.yaml envs/env_scvi.yaml \
        containers/Dockerfile.scanpy containers/Dockerfile.scvi containers/Dockerfile.r
git commit -m "feat: add conda envs and Dockerfiles for scanpy, scvi, r containers"
```

---

## Chunk 2: Per-sample modules

### Task 3: Per-sample Nextflow modules

**Files:**
- Create: `workflow/modules/load.nf`
- Create: `workflow/modules/ambient.nf`
- Create: `workflow/modules/qc.nf`
- Create: `workflow/modules/filter.nf`
- Create: `workflow/modules/doublets.nf`
- Create: `workflow/modules/preprocess.nf`

- [ ] **Step 1: Write `workflow/modules/load.nf`**

```groovy
// workflow/modules/load.nf
process LOAD {
    label 'scanpy'
    tag   "${sample}"

    input:
    tuple val(sample), path(input_dir)

    output:
    tuple val(sample), path("${sample}_loaded.h5ad")

    script:
    """
    scanpy-workflow load \
        --input    ${input_dir} \
        --sample   ${sample} \
        --output   ${sample}_loaded.h5ad
    """
}
```

- [ ] **Step 2: Write `workflow/modules/ambient.nf`**

Dual-dispatch: SoupX (Rscript) or DecontX (Rscript). Both R; no Python CLI path.

```groovy
// workflow/modules/ambient.nf
process AMBIENT {
    label 'r_env'
    tag   "${sample}"

    input:
    tuple val(sample), path(h5ad)
    path  raw_dir      // CellRanger raw_feature_bc_matrix/ (only used by SoupX)
    val   method       // "soupx" | "decontx"

    output:
    tuple val(sample), path("${sample}_ambient.h5ad")

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    if (method == "soupx") {
        """
        Rscript ${script_dir}/ambient/soupx.R \
            --filtered-h5ad ${h5ad} \
            --raw-dir        ${raw_dir} \
            --output         ${sample}_ambient.h5ad
        """
    } else if (method == "decontx") {
        """
        Rscript ${script_dir}/ambient/decontx.R \
            --input  ${h5ad} \
            --output ${sample}_ambient.h5ad
        """
    } else {
        error "ambient.nf: unknown method '${method}'. Expected: soupx | decontx"
    }
}
```

- [ ] **Step 3: Write `workflow/modules/qc.nf`**

```groovy
// workflow/modules/qc.nf
process QC {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/qc", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_qc.h5ad")
    path  "qc_plots/"

    script:
    """
    scanpy-workflow qc \
        --input  ${h5ad} \
        --output ${sample}_qc.h5ad \
        --plots  qc_plots/
    """
}
```

- [ ] **Step 4: Write `workflow/modules/filter.nf`**

```groovy
// workflow/modules/filter.nf
process FILTER {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/filtered", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_filtered.h5ad")

    script:
    """
    scanpy-workflow filter \
        --input        ${h5ad} \
        --output       ${sample}_filtered.h5ad \
        --min-genes    ${params.filter.min_genes} \
        --max-genes    ${params.filter.max_genes} \
        --min-counts   ${params.filter.min_counts} \
        --max-counts   ${params.filter.max_counts} \
        --max-pct-mito ${params.filter.max_pct_mito} \
        --max-pct-ribo ${params.filter.max_pct_ribo} \
        --min-cells    ${params.filter.min_cells}
    """
}
```

- [ ] **Step 5: Write `workflow/modules/doublets.nf`**

```groovy
// workflow/modules/doublets.nf
process DOUBLETS {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/doublets", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_doublets.h5ad")

    script:
    """
    scanpy-workflow doublets \
        --input  ${h5ad} \
        --output ${sample}_doublets.h5ad
    """
}
```

- [ ] **Step 6: Write `workflow/modules/preprocess.nf`**

Two separate processes — one per environment. `per_sample.nf` branches on `params.normalization.method`. Nextflow labels are static and cannot be changed at runtime based on input values.

```groovy
// workflow/modules/preprocess.nf

// library_size normalization: Python CLI, env_scanpy
process NORMALIZE {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/preprocessed", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)
    val   target_sum

    output:
    tuple val(sample), path("${sample}_normalized.h5ad")

    script:
    """
    scanpy-workflow normalize \
        --input      ${h5ad} \
        --output     ${sample}_normalized.h5ad \
        --target-sum ${target_sum}
    """
}

// scran normalization: R script, env_r
process NORMALIZE_SCRAN {
    label     'r_env'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/preprocessed", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_normalized.h5ad")

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/preprocessing/scran.R \
        --input  ${h5ad} \
        --output ${sample}_normalized.h5ad
    """
}

process HVG_PER_SAMPLE {
    label 'scanpy'
    tag   "${sample}"

    input:
    tuple val(sample), path(h5ad)
    val   n_top_genes
    val   flavor

    output:
    tuple val(sample), path("${sample}_hvg.h5ad")

    script:
    """
    scanpy-workflow hvg \
        --input       ${h5ad} \
        --output      ${sample}_hvg.h5ad \
        --mode        per-sample \
        --n-top-genes ${n_top_genes} \
        --flavor      ${flavor}
    """
}
```

- [ ] **Step 7: Commit**

```bash
git add workflow/modules/load.nf workflow/modules/ambient.nf \
        workflow/modules/qc.nf workflow/modules/filter.nf \
        workflow/modules/doublets.nf workflow/modules/preprocess.nf
git commit -m "feat: add per-sample Nextflow modules (load, ambient, qc, filter, doublets, normalize, hvg)"
```

---

### Task 4: Per-sample subworkflow

**Files:**
- Create: `workflow/subworkflows/per_sample.nf`

- [ ] **Step 1: Write `workflow/subworkflows/per_sample.nf`**

```groovy
// workflow/subworkflows/per_sample.nf
include { LOAD          } from '../modules/load'
include { AMBIENT       } from '../modules/ambient'
include { QC            } from '../modules/qc'
include { FILTER        } from '../modules/filter'
include { DOUBLETS      } from '../modules/doublets'
include { NORMALIZE       } from '../modules/preprocess'
include { NORMALIZE_SCRAN  } from '../modules/preprocess'
include { HVG_PER_SAMPLE  } from '../modules/preprocess'

workflow PER_SAMPLE {
    take:
    // Channel of [ sample_name, cellranger_dir ]
    samples_ch  // tuple(val, path)

    main:
    // Step 1: Load Cell Ranger output → .h5ad
    loaded_ch = LOAD(samples_ch)

    // Step 2 (optional): Ambient RNA removal
    if (params.ambient.enabled) {
        // SoupX also needs the raw_feature_bc_matrix dir
        raw_dirs_ch = samples_ch.map { sample, dir ->
            tuple(sample, file("${dir}/raw_feature_bc_matrix"))
        }
        // Join on sample to pair loaded h5ad with its raw dir
        ambient_input_ch = loaded_ch.join(raw_dirs_ch).map { sample, h5ad, raw ->
            tuple(sample, h5ad, raw)
        }
        post_ambient_ch = AMBIENT(
            ambient_input_ch.map { s, h, r -> tuple(s, h) },
            ambient_input_ch.map { s, h, r -> r },
            params.ambient.method
        )
    } else {
        post_ambient_ch = loaded_ch
    }

    // Step 3: QC metrics
    qc_ch = QC(post_ambient_ch)

    // Step 4: Filter
    filtered_ch = FILTER(qc_ch[0])   // qc_ch emits [tuple, path(plots)]

    // Step 5 (optional): Doublet detection
    if (params.doublet.enabled) {
        doublet_ch = DOUBLETS(filtered_ch)
    } else {
        doublet_ch = filtered_ch
    }

    // Step 6: Normalization — branch on method to select correct env
    if (params.normalization.method == "scran") {
        norm_ch = NORMALIZE_SCRAN(doublet_ch)
    } else {
        norm_ch = NORMALIZE(doublet_ch, params.normalization.target_sum)
    }

    // Step 7: Per-sample HVG
    hvg_ch = HVG_PER_SAMPLE(
        norm_ch,
        params.hvg.n_top_genes,
        params.hvg.flavor
    )

    emit:
    preprocessed = hvg_ch   // tuple(sample_name, preprocessed.h5ad)
}
```

- [ ] **Step 2: Commit**

```bash
git add workflow/subworkflows/per_sample.nf
git commit -m "feat: add per_sample subworkflow"
```

---

## Chunk 3: Integration modules + subworkflow

### Task 5: Integration modules

**Files:**
- Create: `workflow/modules/merge.nf`
- Create: `workflow/modules/hvg_postmerge.nf`
- Create: `workflow/modules/scale.nf`
- Create: `workflow/modules/pca.nf`
- Create: `workflow/modules/batch_correction.nf`
- Create: `workflow/modules/imputation.nf`

- [ ] **Step 1: Write `workflow/modules/merge.nf`**

```groovy
// workflow/modules/merge.nf
process MERGE {
    label     'scanpy'
    publishDir "${params.output_dir}/merged", mode: 'copy'

    input:
    path h5ad_files     // collected list of per-sample .h5ad
    val  sample_names   // space-separated string of sample names
    val  batch_key      // obs column for batch identity

    output:
    path "merged.h5ad"

    script:
    """
    scanpy-workflow merge \
        --inputs    ${h5ad_files} \
        --samples   ${sample_names} \
        --batch-key ${batch_key} \
        --output    merged.h5ad
    """
}
```

- [ ] **Step 2: Write `workflow/modules/hvg_postmerge.nf`**

```groovy
// workflow/modules/hvg_postmerge.nf
process HVG_POST_MERGE {
    label     'scanpy'
    publishDir "${params.output_dir}/merged", mode: 'copy'

    input:
    path h5ad
    val  n_top_genes
    val  flavor
    val  batch_key     // passed to batch_key param for scib-style HVG

    output:
    path "merged_hvg.h5ad"

    script:
    """
    scanpy-workflow hvg \
        --input        ${h5ad} \
        --output       merged_hvg.h5ad \
        --mode         post-merge \
        --n-top-genes  ${n_top_genes} \
        --flavor       ${flavor} \
        --batch-key    ${batch_key}
    """
}
```

- [ ] **Step 3: Write `workflow/modules/scale.nf`**

```groovy
// workflow/modules/scale.nf
process SCALE {
    label     'scanpy'
    publishDir "${params.output_dir}/scaled", mode: 'copy'

    input:
    path h5ad

    output:
    path "scaled.h5ad"

    script:
    """
    scanpy-workflow scale \
        --input  ${h5ad} \
        --output scaled.h5ad
    """
}
```

- [ ] **Step 4: Write `workflow/modules/pca.nf`**

```groovy
// workflow/modules/pca.nf
process PCA {
    label     'scanpy'
    publishDir "${params.output_dir}/pca", mode: 'copy'

    input:
    path h5ad
    val  n_comps
    val  scvi_path   // "true" | "false" — skips scale, asserts adata.X == log_norm

    output:
    path "pca.h5ad"

    script:
    """
    scanpy-workflow pca \
        --input     ${h5ad} \
        --output    pca.h5ad \
        --n-comps   ${n_comps} \
        --scvi-path ${scvi_path}
    """
}
```

- [ ] **Step 5: Write `workflow/modules/batch_correction.nf`**

Dual-dispatch: Harmony and BBKNN use `env_scanpy`; scVI uses `env_scvi`.

Two separate processes — one per environment. `integration.nf` branches on `params.batch_correction.method`. Nextflow labels are static.

```groovy
// workflow/modules/batch_correction.nf

// Harmony and BBKNN: Python CLI, env_scanpy
process BATCH_CORRECT {
    label     'scanpy'
    publishDir "${params.output_dir}/batch_correction", mode: 'copy'

    input:
    path h5ad
    val  method      // "harmony" | "bbknn"
    val  batch_key

    output:
    path "batch_corrected.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow batch-correct \
        --input     ${h5ad} \
        --output    batch_corrected.h5ad \
        --method    ${method} \
        --batch-key ${batch_key} \
        --eval-dir  evaluation/
    """
}

// scVI batch correction: Python CLI, env_scvi (requires torch + scvi-tools)
process BATCH_CORRECT_SCVI {
    label     'scvi'
    publishDir "${params.output_dir}/batch_correction", mode: 'copy'

    input:
    path h5ad
    val  batch_key

    output:
    path "batch_corrected.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow batch-correct \
        --input     ${h5ad} \
        --output    batch_corrected.h5ad \
        --method    scvi \
        --batch-key ${batch_key} \
        --eval-dir  evaluation/
    """
}
```

- [ ] **Step 6: Write `workflow/modules/imputation.nf`**

Three separate processes — one per environment. `integration.nf` branches on `params.imputation.method`. Nextflow labels are static.

```groovy
// workflow/modules/imputation.nf

// MAGIC: Python CLI, env_scanpy
process IMPUTE {
    label     'scanpy'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad
    val  evaluate     // "true" | "false"

    output:
    path "imputed.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow impute \
        --input    ${h5ad} \
        --output   imputed.h5ad \
        --method   magic \
        --evaluate ${evaluate} \
        --eval-dir evaluation/
    """
}

// scVI imputation: Python CLI, env_scvi (requires torch + scvi-tools)
process IMPUTE_SCVI {
    label     'scvi'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad
    val  evaluate

    output:
    path "imputed.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow impute \
        --input    ${h5ad} \
        --output   imputed.h5ad \
        --method   scvi \
        --evaluate ${evaluate} \
        --eval-dir evaluation/
    """
}

// ALRA imputation: R script, env_r
process IMPUTE_ALRA {
    label     'r_env'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad

    output:
    path "imputed.h5ad"

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/imputation/alra.R \
        --input  ${h5ad} \
        --output imputed.h5ad
    """
}
```

- [ ] **Step 7: Commit**

```bash
git add workflow/modules/merge.nf workflow/modules/hvg_postmerge.nf \
        workflow/modules/scale.nf workflow/modules/pca.nf \
        workflow/modules/batch_correction.nf workflow/modules/imputation.nf
git commit -m "feat: add integration Nextflow modules (merge, hvg, scale, pca, BC, imputation)"
```

---

### Task 6: Integration subworkflow

**Files:**
- Create: `workflow/subworkflows/integration.nf`

- [ ] **Step 1: Write `workflow/subworkflows/integration.nf`**

```groovy
// workflow/subworkflows/integration.nf
include { MERGE              } from '../modules/merge'
include { HVG_POST_MERGE     } from '../modules/hvg_postmerge'
include { SCALE              } from '../modules/scale'
include { PCA                } from '../modules/pca'
include { BATCH_CORRECT      } from '../modules/batch_correction'
include { BATCH_CORRECT_SCVI } from '../modules/batch_correction'
include { IMPUTE             } from '../modules/imputation'
include { IMPUTE_SCVI        } from '../modules/imputation'
include { IMPUTE_ALRA        } from '../modules/imputation'

workflow INTEGRATION {
    take:
    preprocessed_ch   // Channel of tuple(sample_name, preprocessed.h5ad)

    main:
    // Step 8: Merge samples
    // Use multiMap to split channel once, avoiding double-consumption.
    preprocessed_ch
        .multiMap { s, h ->
            h5ads:   h
            samples: s
        }
        .set { split_ch }

    all_h5ads_ch    = split_ch.h5ads.collect()
    all_samples_str = split_ch.samples.collect().map { names -> names.join(' ') }

    merged_h5ad = MERGE(
        all_h5ads_ch,
        all_samples_str,
        params.batch_correction.batch_key
    )

    // Step 9: Post-merge HVG
    hvg_h5ad = HVG_POST_MERGE(
        merged_h5ad,
        params.hvg.post_merge_n_top_genes,
        params.hvg.flavor,
        params.batch_correction.batch_key
    )

    // Step 10: Scale (skip if scVI batch correction — scVI path keeps adata.X = log_norm)
    boolean scvi_bc = (params.batch_correction.enabled &&
                       params.batch_correction.method == "scvi")
    if (scvi_bc) {
        pre_pca_h5ad = hvg_h5ad
    } else {
        pre_pca_h5ad = SCALE(hvg_h5ad)
    }

    // Step 11: PCA
    pca_h5ad = PCA(pre_pca_h5ad, params.pca.n_comps, scvi_bc.toString())

    // Step 12 (optional): Batch correction — branch on method to select env
    if (params.batch_correction.enabled) {
        if (params.batch_correction.method == "scvi") {
            bc_result    = BATCH_CORRECT_SCVI(pca_h5ad, params.batch_correction.batch_key)
        } else {
            bc_result    = BATCH_CORRECT(pca_h5ad, params.batch_correction.method,
                                         params.batch_correction.batch_key)
        }
        post_bc_h5ad = bc_result[0]
    } else {
        post_bc_h5ad = pca_h5ad
    }

    // Step 13 (optional): Imputation — branch on method to select env
    if (params.imputation.enabled) {
        def eval_str = params.imputation.evaluate.toString()
        if (params.imputation.method == "scvi") {
            imp_result = IMPUTE_SCVI(post_bc_h5ad, eval_str)
        } else if (params.imputation.method == "alra") {
            imp_result = IMPUTE_ALRA(post_bc_h5ad)
        } else {
            imp_result = IMPUTE(post_bc_h5ad, eval_str)
        }
        final_h5ad = imp_result[0]
    } else {
        final_h5ad = post_bc_h5ad
    }

    emit:
    integrated = final_h5ad   // path to integrated .h5ad
}
```

- [ ] **Step 2: Commit**

```bash
git add workflow/subworkflows/integration.nf
git commit -m "feat: add integration subworkflow (merge→HVG→scale→PCA→BC→imputation)"
```

---

## Chunk 4: Downstream modules + main.nf

### Task 7: Downstream modules

**Files:**
- Create: `workflow/modules/clustering.nf`
- Create: `workflow/modules/annotation.nf`
- Create: `workflow/modules/de.nf`
- Create: `workflow/modules/report.nf`

- [ ] **Step 1: Write `workflow/modules/clustering.nf`**

```groovy
// workflow/modules/clustering.nf
process CLUSTER {
    label     'scanpy'
    publishDir "${params.output_dir}/clustering", mode: 'copy'

    input:
    path h5ad
    val  algorithm   // "leiden" | "louvain"
    val  resolution
    val  n_neighbors
    val  evaluate    // "true" | "false"

    output:
    path "clustered.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow cluster \
        --input       ${h5ad} \
        --output      clustered.h5ad \
        --algorithm   ${algorithm} \
        --resolution  ${resolution} \
        --n-neighbors ${n_neighbors} \
        --evaluate    ${evaluate} \
        --eval-dir    evaluation/
    """
}
```

- [ ] **Step 2: Write `workflow/modules/annotation.nf`**

Dispatches rank_genes (always Python) + automated annotation (CellTypist=Python, scType=Rscript) in parallel.

```groovy
// workflow/modules/annotation.nf
process RANK_GENES {
    label     'scanpy'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path h5ad
    val  groupby

    output:
    path "rank_genes.h5ad"
    path "marker_plots/"

    script:
    """
    scanpy-workflow rank-genes \
        --input   ${h5ad} \
        --output  rank_genes.h5ad \
        --groupby ${groupby} \
        --plots   marker_plots/
    """
}

process ANNOTATE_CELLTYPIST {
    label     'scanpy'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path  h5ad
    val   model
    val   groupby

    output:
    path "annotated_celltypist.h5ad"

    script:
    """
    scanpy-workflow annotate \
        --input   ${h5ad} \
        --output  annotated_celltypist.h5ad \
        --model   ${model} \
        --groupby ${groupby}
    """
}

process ANNOTATE_SCTYPE {
    label     'r_env'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path h5ad
    path markers_json
    val  groupby

    output:
    path "annotated_sctype.h5ad"

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/annotation/sctype.R \
        --input   ${h5ad} \
        --markers ${markers_json} \
        --groupby ${groupby} \
        --output  annotated_sctype.h5ad
    """
}
```

- [ ] **Step 3: Write `workflow/modules/de.nf`**

```groovy
// workflow/modules/de.nf
process DE {
    label     'scanpy'
    publishDir "${params.output_dir}/de", mode: 'copy'

    input:
    path h5ad
    val  method    // "wilcoxon" | "t-test" | "logreg"
    val  groupby

    output:
    path "de_results.csv"
    path "de_plots/"

    script:
    """
    scanpy-workflow de \
        --input   ${h5ad} \
        --output  de_results.csv \
        --method  ${method} \
        --groupby ${groupby} \
        --plots   de_plots/
    """
}
```

- [ ] **Step 4: Write `workflow/modules/report.nf`**

```groovy
// workflow/modules/report.nf
process REPORT {
    label     'scanpy'
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
    path h5ad
    path de_csv
    val  output_dir

    output:
    path "report.html"
    path "final.h5ad"

    script:
    """
    scanpy-workflow report \
        --input      ${h5ad} \
        --de-csv     ${de_csv} \
        --output-dir . \
        --output     report.html \
        --final-h5ad final.h5ad
    """
}
```

- [ ] **Step 5: Commit**

```bash
git add workflow/modules/clustering.nf workflow/modules/annotation.nf \
        workflow/modules/de.nf workflow/modules/report.nf
git commit -m "feat: add downstream Nextflow modules (cluster, annotate, de, report)"
```

---

### Task 8: main.nf

**Files:**
- Create: `workflow/main.nf`

- [ ] **Step 1: Write `workflow/main.nf`**

```groovy
#!/usr/bin/env nextflow
// workflow/main.nf — Scanpy end-to-end single-cell workflow
// DSL2 is the default in Nextflow >= 22.03

include { PER_SAMPLE  } from './subworkflows/per_sample'
include { INTEGRATION } from './subworkflows/integration'
include { CLUSTER     } from './modules/clustering'
include { RANK_GENES  } from './modules/annotation'
include { ANNOTATE_CELLTYPIST } from './modules/annotation'
include { ANNOTATE_SCTYPE     } from './modules/annotation'
include { DE          } from './modules/de'
include { REPORT      } from './modules/report'

// Note: workflow/lib/validate.groovy is auto-loaded by Nextflow; no include needed.

workflow {
    // Config validation (exits with error on invalid combos)
    // validateConfig() is defined in workflow/lib/validate.groovy (auto-loaded)
    validateConfig(params)

    // Build per-sample channel: tuple(sample_name, cellranger_dir)
    samples_ch = Channel.fromList(params.samples)
        .map { sample -> tuple(sample, file("${params.input_dir}/${sample}")) }

    // Steps 1-7: Per-sample processing (parallelized)
    PER_SAMPLE(samples_ch)

    // Steps 8-13: Integration (merge → HVG → scale → PCA → BC → imputation)
    INTEGRATION(PER_SAMPLE.out.preprocessed)

    // Steps 14-15: Neighbors + clustering + evaluation
    cluster_result = CLUSTER(
        INTEGRATION.out.integrated,
        params.clustering.algorithm,
        params.clustering.resolution,
        params.clustering.n_neighbors,
        params.clustering.evaluate.toString()
    )
    // Convert to value channel so multiple downstream processes can consume it.
    // Queue channels are single-consumer; .first() re-emits to each subscriber.
    clustered_h5ad = cluster_result[0].first()

    // Step 16a: Unsupervised annotation (rank genes)
    if (params.annotation.unsupervised) {
        rank_result = RANK_GENES(
            clustered_h5ad,
            params.de.groupby ?: params.clustering.algorithm
        )
        annotated_h5ad = rank_result[0]
    } else {
        annotated_h5ad = clustered_h5ad
    }

    // Step 16b: Automated annotation (CellTypist or scType) — runs in parallel with 16a
    // validateConfig() already guarantees sctype_markers is set if method == "sctype"
    if (params.annotation.automated) {
        if (params.annotation.method == "celltypist") {
            ANNOTATE_CELLTYPIST(
                clustered_h5ad,
                params.annotation.celltypist_model,
                params.clustering.algorithm
            )
        } else if (params.annotation.method == "sctype") {
            ANNOTATE_SCTYPE(
                clustered_h5ad,
                file(params.annotation.sctype_markers),
                params.clustering.algorithm
            )
        }
    }

    // Step 17: Differential expression
    de_result = DE(
        annotated_h5ad,
        params.de.method,
        params.de.groupby ?: params.clustering.algorithm
    )

    // Step 18: HTML report + final .h5ad
    if (params.report) {
        REPORT(
            annotated_h5ad,
            de_result[0],
            params.output_dir
        )
    }
}
```

- [ ] **Step 2: Write `workflow/lib/validate.groovy`**

```groovy
// workflow/lib/validate.groovy
// Config validation called from main.nf at workflow startup.

def validateConfig(params) {
    // Hard error: scVI cannot be used for both batch correction and imputation
    boolean bc_scvi  = params.batch_correction.enabled &&
                       params.batch_correction.method == "scvi"
    boolean imp_scvi = params.imputation.enabled &&
                       params.imputation.method == "scvi"
    if (bc_scvi && imp_scvi) {
        error """
        Invalid configuration: scVI cannot be used for both batch correction and imputation simultaneously.
        Set batch_correction.method to 'harmony' or 'bbknn', or set imputation.method to 'magic' or 'alra'.
        """
    }

    // Hard error: scType requires a markers file
    if (params.annotation.automated &&
        params.annotation.method == "sctype" &&
        params.annotation.sctype_markers == null) {
        error "Invalid configuration: annotation.sctype_markers must be set when annotation.method is 'sctype'."
    }

    // Warning: de.groupby column won't exist at startup (created by clustering)
    if (params.de.groupby && params.de.groupby != params.clustering.algorithm) {
        log.warn "de.groupby ('${params.de.groupby}') differs from clustering.algorithm ('${params.clustering.algorithm}'). Ensure the column will be created before DE step."
    }
}
```

- [ ] **Step 3: Commit**

```bash
git add workflow/main.nf workflow/lib/validate.groovy
git commit -m "feat: add main.nf orchestration and config validation"
```

---

## Chunk 5: Integration tests + CI update

### Task 9: Nextflow integration tests

**Files:**
- Create: `tests/integration/test_pipeline_smoke.py`
- Modify: `tests/conftest.py` (add synthetic CellRanger fixture dirs)

- [ ] **Step 1: Write synthetic CellRanger fixture to `tests/conftest.py`**

The integration smoke test needs a minimal CellRanger-like directory for 2 samples. Add to `tests/conftest.py`:

```python
# tests/conftest.py  (append to existing conftest)
import gzip
import io
import os
import scipy.io
import scipy.sparse
import numpy as np
import pytest

def write_cellranger_dir(base_dir: str, n_cells: int = 50, n_genes: int = 100,
                          sample: str = "sample") -> str:
    """
    Write a minimal CellRanger output structure:
      base_dir/<sample>/filtered_feature_bc_matrix/{matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz}
      base_dir/<sample>/raw_feature_bc_matrix/     {same files, more barcodes}
    """
    rng = np.random.default_rng(42)
    barcodes = [f"ACGT{i:04d}-1" for i in range(n_cells)]
    genes    = [f"GENE{i:04d}" for i in range(n_genes)]
    counts   = rng.poisson(2, size=(n_genes, n_cells)).astype("float32")
    mat      = scipy.sparse.csc_matrix(counts)

    for subdir, bc_list in [
        ("filtered_feature_bc_matrix", barcodes),
        ("raw_feature_bc_matrix",      barcodes + [f"TTTT{i:04d}-1" for i in range(200)])
    ]:
        d = os.path.join(base_dir, sample, subdir)
        os.makedirs(d, exist_ok=True)

        # matrix.mtx.gz
        buf = io.BytesIO()
        scipy.io.mmwrite(buf, mat if subdir.startswith("filtered") else
                         scipy.sparse.csc_matrix(rng.poisson(0.1, size=(n_genes, len(bc_list)))))
        with gzip.open(os.path.join(d, "matrix.mtx.gz"), "wb") as f:
            f.write(buf.getvalue())

        # barcodes.tsv.gz
        with gzip.open(os.path.join(d, "barcodes.tsv.gz"), "wt") as f:
            f.write("\n".join(bc_list) + "\n")

        # features.tsv.gz
        with gzip.open(os.path.join(d, "features.tsv.gz"), "wt") as f:
            for g in genes:
                f.write(f"{g}\t{g}\tGene Expression\n")

    return os.path.join(base_dir, sample)


@pytest.fixture(scope="session")
def cellranger_dir(tmp_path_factory):
    """Two-sample synthetic CellRanger fixture for integration tests."""
    base = str(tmp_path_factory.mktemp("cellranger"))
    write_cellranger_dir(base, sample="sample_A")
    write_cellranger_dir(base, sample="sample_B")
    return base
```

- [ ] **Step 2: Write `tests/integration/test_pipeline_smoke.py`**

```python
# tests/integration/test_pipeline_smoke.py
"""
Smoke test: run the Nextflow pipeline on a synthetic 2-sample dataset.
Verifies:
  - Pipeline exits 0
  - Key output files exist (final.h5ad, de_results.csv, report.html)
  - final.h5ad has the expected AnnData structure
Skipped if Nextflow is not on PATH (CI installs it; local dev may skip).
"""
import os
import subprocess
import shutil
import pytest
import anndata as ad


def nextflow_available():
    return shutil.which("nextflow") is not None


@pytest.mark.skipif(not nextflow_available(), reason="Nextflow not on PATH")
def test_pipeline_smoke(cellranger_dir, tmp_path):
    out_dir = str(tmp_path / "results")
    params  = str(tmp_path / "params.yaml")

    # Write minimal params.yaml for smoke run
    with open(params, "w") as f:
        f.write(f"""\
input_dir: "{cellranger_dir}"
samples: [sample_A, sample_B]
output_dir: "{out_dir}"
ambient:
  enabled: false
doublet:
  enabled: false
normalization:
  method: library_size
  target_sum: 10000
hvg:
  n_top_genes: 50
  post_merge_n_top_genes: 30
  flavor: seurat_v3
pca:
  n_comps: 10
batch_correction:
  enabled: true
  method: harmony
  batch_key: sample
imputation:
  enabled: false
clustering:
  algorithm: leiden
  resolution: 0.5
  n_neighbors: 5
  evaluate: false
annotation:
  unsupervised: true
  automated: false
de:
  method: wilcoxon
  groupby: leiden
report: true
""")

    # Resolve repo root so relative paths in main.nf work regardless of pytest invocation dir
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    result = subprocess.run(
        ["nextflow", "run", "workflow/main.nf",
         "-params-file", params,
         "-profile", "conda",
         "-resume"],
        capture_output=True, text=True, timeout=600,
        cwd=repo_root
    )

    assert result.returncode == 0, (
        f"Nextflow pipeline failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Key output files
    assert os.path.exists(os.path.join(out_dir, "report", "final.h5ad")), \
        "Missing final.h5ad"
    assert os.path.exists(os.path.join(out_dir, "de", "de_results.csv")), \
        "Missing de_results.csv"
    assert os.path.exists(os.path.join(out_dir, "report", "report.html")), \
        "Missing report.html"

    # Structural checks on final.h5ad
    adata = ad.read_h5ad(os.path.join(out_dir, "report", "final.h5ad"))
    assert "leiden" in adata.obs.columns, "Missing leiden column in obs"
    assert "counts" in adata.layers,      "Missing counts layer"
    assert "log_norm" in adata.layers,    "Missing log_norm layer"

    # DE CSV has required columns
    import pandas as pd
    de = pd.read_csv(os.path.join(out_dir, "de", "de_results.csv"))
    required_cols = {"group", "gene", "score", "logfoldchange", "pval", "pval_adj"}
    assert required_cols.issubset(set(de.columns)), \
        f"DE CSV missing columns: {required_cols - set(de.columns)}"
```

- [ ] **Step 3: Commit**

```bash
git add tests/conftest.py tests/integration/test_pipeline_smoke.py
git commit -m "test: add integration smoke test with synthetic CellRanger fixture"
```

---

### Task 10: Update CI for Nextflow integration tests

**Files:**
- Modify: `.github/workflows/ci.yml` (already created in Plan 2)

- [ ] **Step 1: Update the `integration` job in `.github/workflows/ci.yml`**

Replace the integration job body to install Nextflow:

```yaml
  integration:
    name: Integration tests (Nextflow smoke)
    runs-on: ubuntu-latest
    if: github.event_name == 'pull_request'
    needs: [python-unit, r-unit]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - uses: nf-core/setup-nextflow@v2
        with:
          version: "23.10.0"
      - uses: conda-incubator/setup-miniconda@v3
        with:
          auto-activate-base: false
      - name: Create scanpy conda env
        run: conda env create -f envs/env_scanpy.yaml
      - name: Install Python package
        run: conda run -n env_scanpy pip install -e ".[dev]"
      - name: Run integration tests
        run: |
          conda run -n env_scanpy pytest tests/integration/ -v --tb=short -x
```

- [ ] **Step 2: Verify YAML is valid**

```bash
python -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))" && echo "YAML valid"
```

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: update integration job to install Nextflow and run smoke test"
```

---

## Final: Verify full pipeline structure

- [ ] **Step 1: Verify workflow directory structure**

```bash
find workflow/ -name "*.nf" | sort
```

Expected output:
```
workflow/main.nf
workflow/modules/ambient.nf
workflow/modules/annotation.nf
workflow/modules/batch_correction.nf
workflow/modules/clustering.nf
workflow/modules/de.nf
workflow/modules/doublets.nf
workflow/modules/filter.nf
workflow/modules/hvg_postmerge.nf
workflow/modules/imputation.nf
workflow/modules/load.nf
workflow/modules/merge.nf
workflow/modules/pca.nf
workflow/modules/preprocess.nf
workflow/modules/report.nf
workflow/modules/scale.nf
workflow/subworkflows/integration.nf
workflow/subworkflows/per_sample.nf
```

- [ ] **Step 2: Syntax check Nextflow files**

```bash
nextflow inspect workflow/main.nf -params-file params.yaml 2>&1 | head -20
```

Expected: No errors; lists workflow structure.

- [ ] **Step 3: Final commit**

```bash
git add -A
git commit -m "feat: complete Plan 3 — Nextflow DSL2 pipeline scaffold"
```
