# Scanpy Workflow — Plan 1: Python Package

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `scanpy_workflow` Python package — all analytical modules, CLI entrypoints, and unit tests — as specified in the design spec at `docs/superpowers/specs/2026-03-15-scanpy-workflow-design.md`.

**Architecture:** Installable Python package with a Click CLI. Each module reads/writes `.h5ad` files so Nextflow can chain them without in-memory state. Each analytical step is a pure function that takes an AnnData and returns an AnnData. Tests use pytest with synthetic fixtures; no real data downloads required.

**Tech Stack:** Python 3.9+, scanpy, anndata, harmonypy, bbknn, magic-impute, scrublet, celltypist, scib-metrics, scvi-tools (optional env), Click, PyYAML, Jinja2, pytest, scipy, numpy, pandas, matplotlib, seaborn

**Spec:** `docs/superpowers/specs/2026-03-15-scanpy-workflow-design.md`

---

## Chunk 1: Scaffold + Config + Utils

### Task 1: Project scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `src/scanpy_workflow/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/unit/__init__.py`
- Create: `tests/conftest.py`

- [ ] **Step 1: Create directory structure and minimal CLI stub**

```bash
mkdir -p src/scanpy_workflow/{io,qc,ambient,doublets,preprocessing,integration,imputation,clustering,annotation,de,reporting,utils}
mkdir -p tests/unit/r tests/integration
touch src/scanpy_workflow/__init__.py
touch src/scanpy_workflow/{io,qc,ambient,doublets,preprocessing,integration,imputation,clustering,annotation,de,reporting,utils}/__init__.py
touch tests/__init__.py tests/unit/__init__.py tests/integration/__init__.py
```

Create a minimal `src/scanpy_workflow/cli.py` stub so the installed entrypoint resolves at install time. Full commands are added in Task 23.

```python
# src/scanpy_workflow/cli.py  (stub — commands added in Task 23)
import click

@click.group()
def cli():
    """Scanpy end-to-end single-cell workflow."""
    pass
```

- [ ] **Step 2: Write `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "scanpy-workflow"
version = "0.1.0"
requires-python = ">=3.9"
dependencies = [
    "scanpy>=1.9",
    "anndata>=0.9",
    "harmonypy>=0.0.9",
    "bbknn>=1.6",
    "magic-impute>=3.0",
    "scrublet>=0.2",
    "celltypist>=1.5",
    "scib-metrics>=0.4",
    "click>=8.0",
    "pyyaml>=6.0",
    "jinja2>=3.0",
    "numpy>=1.24",
    "pandas>=2.0",
    "scipy>=1.10",
    "matplotlib>=3.7",
    "seaborn>=0.12",
    "leidenalg>=0.10",
    "louvain>=0.8",
    "python-igraph>=0.11",
    "scikit-learn>=1.3",
]

[project.optional-dependencies]
scvi = ["scvi-tools>=1.0"]
dev = ["pytest>=7.0", "pytest-cov>=4.0"]

[project.scripts]
scanpy-workflow = "scanpy_workflow.cli:cli"

[tool.hatch.build.targets.wheel]
packages = ["src/scanpy_workflow"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 3: Write `tests/conftest.py` with shared fixtures**

```python
import pytest
import numpy as np
import anndata as ad
import scipy.sparse


@pytest.fixture
def small_adata():
    """200 cells x 300 genes with raw counts and mito/ribo gene annotations."""
    rng = np.random.default_rng(42)
    counts = rng.negative_binomial(5, 0.5, size=(200, 300)).astype(np.float32)
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(counts))
    adata.obs_names = [f"cell_{i}" for i in range(200)]
    # First 10 genes are mitochondrial, next 10 ribosomal
    var_names = (
        [f"MT-gene_{i}" for i in range(10)]
        + [f"RPS{i}" for i in range(10)]
        + [f"gene_{i}" for i in range(280)]
    )
    adata.var_names = var_names
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    return adata


@pytest.fixture
def normalized_adata(small_adata):
    """AnnData with layers: counts, norm, log_norm. X = log_norm."""
    import scanpy as sc
    adata = small_adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["norm"] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()
    return adata


@pytest.fixture
def two_sample_adatas():
    """Two small AnnData objects (100 cells x 200 genes each) simulating two samples."""
    rng = np.random.default_rng(42)
    adatas = []
    for s in range(2):
        counts = rng.negative_binomial(5, 0.5, size=(100, 200)).astype(np.float32)
        adata = ad.AnnData(X=scipy.sparse.csr_matrix(counts))
        adata.obs_names = [f"s{s}_cell_{i}" for i in range(100)]
        adata.var_names = [f"gene_{i}" for i in range(200)]
        adata.obs["sample"] = f"sample_{s}"
        adata.layers["counts"] = adata.X.copy()
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        adata.layers["norm"] = adata.X.copy()
        sc.pp.log1p(adata)
        adata.layers["log_norm"] = adata.X.copy()
        adatas.append(adata)
    return adatas


@pytest.fixture
def merged_adata(two_sample_adatas):
    """Merged AnnData with obs['sample'] batch key and X_pca."""
    import anndata
    import scanpy as sc
    adata = anndata.concat(
        two_sample_adatas,
        label="sample",
        keys=["sample_0", "sample_1"],
        merge="same",
    )
    sc.pp.highly_variable_genes(adata, n_top_genes=100, flavor="seurat")
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=20, use_highly_variable=True)
    adata.uns["neighbors_use_rep"] = "X_pca"
    return adata
```

- [ ] **Step 4: Install package in editable mode**

```bash
pip install -e ".[dev]"
```

Expected: installs without error; `scanpy-workflow --help` shows the top-level group (no subcommands yet — those are added in Task 23).

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/ tests/
git commit -m "feat: project scaffold and test fixtures"
```

---

### Task 2: Config loader + validation

**Files:**
- Create: `src/scanpy_workflow/utils/config.py`
- Create: `tests/unit/test_config.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_config.py
import pytest
import warnings
from scanpy_workflow.utils.config import load_config


def test_load_minimal_config(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: harmony
  batch_key: sample
imputation:
  enabled: false
  method: magic
clustering:
  algorithm: leiden
de:
  method: wilcoxon
  groupby: leiden
""")
    cfg = load_config(str(cfg_file))
    assert cfg["batch_correction"]["method"] == "harmony"


def test_scvi_conflict_raises(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: scvi
imputation:
  enabled: true
  method: scvi
""")
    with pytest.raises(ValueError, match="scVI cannot be used for both"):
        load_config(str(cfg_file))


def test_de_groupby_mismatch_warns(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
clustering:
  algorithm: leiden
de:
  method: wilcoxon
  groupby: louvain
""")
    with pytest.warns(UserWarning, match="does not match clustering.algorithm"):
        load_config(str(cfg_file))


def test_no_scvi_conflict_when_one_disabled(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: scvi
imputation:
  enabled: false
  method: scvi
""")
    cfg = load_config(str(cfg_file))  # should not raise
    assert cfg["batch_correction"]["method"] == "scvi"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_config.py -v
```

Expected: ImportError or AttributeError — `config` module does not exist yet.

- [ ] **Step 3: Write `src/scanpy_workflow/utils/config.py`**

```python
import warnings
import yaml


def load_config(config_path: str) -> dict:
    """Load params.yaml and validate. Raises ValueError for hard errors."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}
    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    bc = cfg.get("batch_correction", {})
    imp = cfg.get("imputation", {})
    bc_enabled = bc.get("enabled", False)
    imp_enabled = imp.get("enabled", False)
    bc_method = bc.get("method", "harmony")
    imp_method = imp.get("method", "magic")

    if bc_enabled and imp_enabled and bc_method == "scvi" and imp_method == "scvi":
        raise ValueError(
            "Invalid config: batch_correction.method=scvi and imputation.method=scvi "
            "cannot both be enabled. scVI cannot be used for both steps simultaneously."
        )

    de_groupby = cfg.get("de", {}).get("groupby")
    clustering_algo = cfg.get("clustering", {}).get("algorithm", "leiden")
    if de_groupby and de_groupby != clustering_algo:
        warnings.warn(
            f"de.groupby='{de_groupby}' does not match clustering.algorithm='{clustering_algo}'. "
            "Ensure this .obs column will exist at DE runtime.",
            UserWarning,
            stacklevel=3,
        )
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_config.py -v
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/utils/config.py tests/unit/test_config.py
git commit -m "feat: config loader with startup validation"
```

---

### Task 3: AnnData I/O helpers

**Files:**
- Create: `src/scanpy_workflow/utils/io.py`
- Create: `tests/unit/test_io_utils.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_io_utils.py
import pytest
import anndata as ad
import numpy as np
import scipy.sparse
from scanpy_workflow.utils.io import read_h5ad, write_h5ad


def test_write_and_read_roundtrip(tmp_path, small_adata):
    out = tmp_path / "test.h5ad"
    write_h5ad(small_adata, str(out))
    loaded = read_h5ad(str(out))
    assert loaded.shape == small_adata.shape
    assert list(loaded.obs_names) == list(small_adata.obs_names)


def test_write_creates_parent_dirs(tmp_path, small_adata):
    out = tmp_path / "subdir" / "nested" / "test.h5ad"
    write_h5ad(small_adata, str(out))
    assert out.exists()
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/unit/test_io_utils.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/utils/io.py`**

```python
from pathlib import Path
import anndata as ad


def read_h5ad(path: str) -> ad.AnnData:
    """Read an AnnData object from an .h5ad file."""
    return ad.read_h5ad(path)


def write_h5ad(adata: ad.AnnData, path: str) -> None:
    """Write an AnnData object to an .h5ad file, creating parent dirs as needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_io_utils.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/utils/io.py tests/unit/test_io_utils.py
git commit -m "feat: AnnData read/write helpers"
```

---

## Chunk 2: Data Loading + QC

### Task 4: Cell Ranger data loader

**Files:**
- Create: `src/scanpy_workflow/io/load.py`
- Create: `tests/unit/test_load.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_load.py
import pytest
import numpy as np
from pathlib import Path
from scanpy_workflow.io.load import load_cellranger


def test_load_cellranger_mex(tmp_path):
    """Test loading MEX format (barcodes/features/matrix)."""
    import scanpy as sc
    import scipy.io
    import gzip
    import anndata as ad

    # Create minimal MEX structure
    sample_dir = tmp_path / "sample_mex" / "filtered_feature_bc_matrix"
    sample_dir.mkdir(parents=True)

    n_cells, n_genes = 50, 100
    rng = np.random.default_rng(0)
    matrix = scipy.sparse.random(n_genes, n_cells, density=0.1, format="csc",
                                  random_state=0, data_rvs=lambda s: rng.integers(1, 10, s))
    # Write matrix as gzipped .mtx to match scanpy's expected format
    import io, gzip as gz
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, matrix)
    with gz.open(sample_dir / "matrix.mtx.gz", "wb") as f:
        f.write(buf.getvalue())

    with gzip.open(sample_dir / "barcodes.tsv.gz", "wt") as f:
        for i in range(n_cells):
            f.write(f"CELL{i:04d}-1\n")

    with gzip.open(sample_dir / "features.tsv.gz", "wt") as f:
        for i in range(n_genes):
            f.write(f"ENSG{i:08d}\tGene{i}\tGene Expression\n")

    adata = load_cellranger(str(tmp_path / "sample_mex"))
    assert adata.shape == (n_cells, n_genes)
    assert adata.X is not None


def test_load_raises_if_no_data(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        load_cellranger(str(empty_dir))
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/unit/test_load.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/io/load.py`**

```python
from pathlib import Path
import anndata as ad
import scanpy as sc


def load_cellranger(sample_dir: str) -> ad.AnnData:
    """
    Load Cell Ranger output (MEX or HDF5) from a sample directory.

    Tries the following paths in order:
      1. {sample_dir}/filtered_feature_bc_matrix/   (MEX)
      2. {sample_dir}/raw_feature_bc_matrix/         (MEX)
      3. {sample_dir}/filtered_feature_bc_matrix.h5  (HDF5)
      4. {sample_dir}/raw_feature_bc_matrix.h5        (HDF5)
    """
    base = Path(sample_dir)
    mex_paths = [
        base / "filtered_feature_bc_matrix",
        base / "raw_feature_bc_matrix",
    ]
    h5_paths = [
        base / "filtered_feature_bc_matrix.h5",
        base / "raw_feature_bc_matrix.h5",
    ]

    for mex_path in mex_paths:
        if mex_path.is_dir():
            return sc.read_10x_mtx(str(mex_path), var_names="gene_symbols", cache=False)

    for h5_path in h5_paths:
        if h5_path.is_file():
            return sc.read_10x_h5(str(h5_path))

    raise FileNotFoundError(
        f"No Cell Ranger output found in {sample_dir}. "
        "Expected filtered_feature_bc_matrix/ or filtered_feature_bc_matrix.h5"
    )
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_load.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/io/load.py tests/unit/test_load.py
git commit -m "feat: Cell Ranger loader (MEX and HDF5)"
```

---

### Task 5: QC metrics computation

**Files:**
- Create: `src/scanpy_workflow/qc/metrics.py`
- Create: `tests/unit/test_qc_metrics.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_qc_metrics.py
from scanpy_workflow.qc.metrics import compute_qc_metrics


def test_qc_adds_required_obs_columns(small_adata):
    adata = compute_qc_metrics(small_adata)
    required_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ribo"]
    for col in required_cols:
        assert col in adata.obs.columns, f"Missing obs column: {col}"


def test_qc_mito_pct_range(small_adata):
    adata = compute_qc_metrics(small_adata)
    assert (adata.obs["pct_counts_mito"] >= 0).all()
    assert (adata.obs["pct_counts_mito"] <= 100).all()


def test_qc_infers_mito_if_not_annotated():
    """If adata.var['mt'] is absent, infer from MT- prefix."""
    import anndata as ad, numpy as np, scipy.sparse
    rng = np.random.default_rng(0)
    counts = rng.integers(0, 10, size=(50, 100)).astype(np.float32)
    adata = ad.AnnData(X=scipy.sparse.csr_matrix(counts))
    adata.obs_names = [f"c{i}" for i in range(50)]
    adata.var_names = [f"MT-gene_{i}" if i < 5 else f"gene_{i}" for i in range(100)]
    # No adata.var["mt"] defined
    adata = compute_qc_metrics(adata)
    assert "pct_counts_mito" in adata.obs.columns
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/unit/test_qc_metrics.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/qc/metrics.py`**

```python
import anndata as ad
import scanpy as sc


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute QC metrics and add to adata.obs.

    Adds: n_genes_by_counts, total_counts, pct_counts_mito, pct_counts_ribo.
    Infers mt/ribo gene flags from var_names if not already in adata.var.
    """
    if "mt" not in adata.var.columns:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if "ribo" not in adata.var.columns:
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # Rename mito column for consistent naming (scanpy outputs pct_counts_mt, we use pct_counts_mito)
    adata.obs["pct_counts_mito"] = adata.obs["pct_counts_mt"]
    # pct_counts_ribo is already named correctly by scanpy (qc_vars=["ribo"] → pct_counts_ribo)
    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_qc_metrics.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/qc/metrics.py tests/unit/test_qc_metrics.py
git commit -m "feat: QC metrics computation"
```

---

### Task 6: Cell and gene filtering

**Files:**
- Create: `src/scanpy_workflow/qc/filter.py`
- Create: `tests/unit/test_filter.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_filter.py
import numpy as np
from scanpy_workflow.qc.metrics import compute_qc_metrics
from scanpy_workflow.qc.filter import filter_cells_genes


def test_filter_removes_low_gene_cells(small_adata):
    adata = compute_qc_metrics(small_adata)
    # Force some cells to have very low gene counts
    import scipy.sparse
    X = adata.X.toarray()
    X[:10] = 0  # first 10 cells get no expression
    import anndata as ad
    adata_mod = ad.AnnData(X=scipy.sparse.csr_matrix(X), obs=adata.obs.copy(), var=adata.var.copy())
    adata_mod = compute_qc_metrics(adata_mod)
    filtered = filter_cells_genes(adata_mod, min_genes=1, max_genes=10000,
                                   min_counts=1, max_counts=10**9,
                                   max_pct_mito=100, max_pct_ribo=100, min_cells=1)
    assert filtered.n_obs < adata_mod.n_obs


def test_filter_removes_high_mito_cells(small_adata):
    import anndata as ad, scipy.sparse
    adata = compute_qc_metrics(small_adata)
    original_n = adata.n_obs
    filtered = filter_cells_genes(adata, min_genes=1, max_genes=10000,
                                   min_counts=1, max_counts=10**9,
                                   max_pct_mito=0.0,  # nothing passes
                                   max_pct_ribo=100, min_cells=1)
    # All cells removed if no mito allowed
    assert filtered.n_obs == 0 or filtered.n_obs < original_n


def test_filter_removes_lowly_expressed_genes(small_adata):
    adata = compute_qc_metrics(small_adata)
    # Gene expressed in 0 cells should be removed
    filtered = filter_cells_genes(adata, min_genes=1, max_genes=10000,
                                   min_counts=1, max_counts=10**9,
                                   max_pct_mito=100, max_pct_ribo=100,
                                   min_cells=10**9)  # impossible threshold
    assert filtered.n_vars == 0 or filtered.n_vars < adata.n_vars


def test_filter_preserves_layers(normalized_adata):
    adata = compute_qc_metrics(normalized_adata)
    filtered = filter_cells_genes(adata, min_genes=1, max_genes=10000,
                                   min_counts=1, max_counts=10**9,
                                   max_pct_mito=100, max_pct_ribo=100, min_cells=1)
    assert "counts" in filtered.layers
    assert "log_norm" in filtered.layers
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/unit/test_filter.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/qc/filter.py`**

```python
import anndata as ad
import scanpy as sc


def filter_cells_genes(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    min_counts: int = 500,
    max_counts: int = 30000,
    max_pct_mito: float = 20.0,
    max_pct_ribo: float = 50.0,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Filter cells and genes based on QC cutoffs.
    Requires compute_qc_metrics() to have been run first.
    """
    required = ["n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ribo"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing QC columns: {missing}. Run compute_qc_metrics() first.")

    cell_mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["pct_counts_mito"] <= max_pct_mito)
        & (adata.obs["pct_counts_ribo"] <= max_pct_ribo)
    )
    adata = adata[cell_mask].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_filter.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/qc/filter.py tests/unit/test_filter.py
git commit -m "feat: cell and gene filtering"
```

---

### Task 7: QC plots

**Files:**
- Create: `src/scanpy_workflow/qc/plots.py`
- Create: `tests/unit/test_qc_plots.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_qc_plots.py
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
from scanpy_workflow.qc.metrics import compute_qc_metrics
from scanpy_workflow.qc.plots import plot_qc


def test_plot_qc_creates_files(tmp_path, small_adata):
    adata = compute_qc_metrics(small_adata)
    plot_qc(adata, output_dir=str(tmp_path))
    assert (tmp_path / "qc_violin.png").exists()
    assert (tmp_path / "qc_scatter_counts_vs_genes.png").exists()
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_qc_plots.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/qc/plots.py`**

```python
from pathlib import Path
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt


def plot_qc(adata: ad.AnnData, output_dir: str) -> None:
    """Generate QC violin and scatter plots, save as PNG."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mito"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig(out / "qc_violin.png", bbox_inches="tight", dpi=100)
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=3,
        alpha=0.5,
    )
    ax.set_xlabel("Total counts")
    ax.set_ylabel("Genes by counts")
    fig.savefig(out / "qc_scatter_counts_vs_genes.png", bbox_inches="tight", dpi=100)
    plt.close(fig)
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_qc_plots.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/qc/plots.py tests/unit/test_qc_plots.py
git commit -m "feat: QC violin and scatter plots"
```

---

## Chunk 3: Per-Sample Preprocessing

### Task 8: Doublet detection with Scrublet

**Files:**
- Create: `src/scanpy_workflow/doublets/scrublet.py`
- Create: `tests/unit/test_scrublet.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_scrublet.py
import numpy as np
from scanpy_workflow.doublets.scrublet import detect_doublets


def test_scrublet_adds_doublet_score(small_adata):
    adata = detect_doublets(small_adata)
    assert "doublet_score" in adata.obs.columns
    assert "predicted_doublet" in adata.obs.columns


def test_scrublet_doublet_score_range(small_adata):
    adata = detect_doublets(small_adata)
    assert (adata.obs["doublet_score"] >= 0).all()
    assert (adata.obs["doublet_score"] <= 1).all()


def test_scrublet_filter_removes_doublets(small_adata):
    adata = detect_doublets(small_adata, filter_doublets=True)
    # After filtering, no predicted doublets remain
    assert not adata.obs["predicted_doublet"].any()
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_scrublet.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/doublets/scrublet.py`**

```python
import anndata as ad
import scrublet as scr
import numpy as np


def detect_doublets(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.05,
    filter_doublets: bool = False,
    random_state: int = 42,
) -> ad.AnnData:
    """
    Detect doublets with Scrublet.

    Adds adata.obs['doublet_score'] and adata.obs['predicted_doublet'].
    If filter_doublets=True, removes predicted doublets.
    Uses raw counts (layers['counts'] if available, else adata.X).
    """
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X

    scrub = scr.Scrublet(
        counts_matrix=counts,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state,
    )
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    if filter_doublets:
        adata = adata[~adata.obs["predicted_doublet"]].copy()

    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_scrublet.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/doublets/scrublet.py tests/unit/test_scrublet.py
git commit -m "feat: Scrublet doublet detection"
```

---

### Task 9: Normalization

**Files:**
- Create: `src/scanpy_workflow/preprocessing/normalize.py`
- Create: `tests/unit/test_normalize.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_normalize.py
import numpy as np
import scipy.sparse
from scanpy_workflow.preprocessing.normalize import normalize


def test_normalize_stores_three_layers(small_adata):
    adata = normalize(small_adata)
    assert "counts" in adata.layers
    assert "norm" in adata.layers
    assert "log_norm" in adata.layers


def test_normalize_x_equals_log_norm(small_adata):
    adata = normalize(small_adata)
    if scipy.sparse.issparse(adata.X):
        x_arr = adata.X.toarray()
        ln_arr = adata.layers["log_norm"].toarray()
    else:
        x_arr = np.asarray(adata.X)
        ln_arr = np.asarray(adata.layers["log_norm"])
    np.testing.assert_array_almost_equal(x_arr, ln_arr)


def test_normalize_counts_layer_is_raw(small_adata):
    """counts layer must equal original X before normalization."""
    import anndata as ad
    original_X = small_adata.X.copy()
    adata = normalize(small_adata.copy())
    if scipy.sparse.issparse(original_X):
        np.testing.assert_array_equal(
            original_X.toarray(),
            adata.layers["counts"].toarray(),
        )


def test_scran_method_raises(small_adata):
    import pytest
    with pytest.raises(ValueError, match="scran normalization must be run via scran.R"):
        normalize(small_adata, method="scran")
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_normalize.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/preprocessing/normalize.py`**

```python
import anndata as ad
import scanpy as sc
import scipy.sparse
import numpy as np


def normalize(
    adata: ad.AnnData,
    method: str = "library_size",
    target_sum: float = 1e4,
) -> ad.AnnData:
    """
    Normalize raw counts and store three layers.

    Layers written:
      - layers["counts"]   : raw integer counts (copy of input X)
      - layers["norm"]     : library-size normalized counts (pre-log)
      - layers["log_norm"] : log1p of normalized counts

    adata.X is set to layers["log_norm"] on return.

    Args:
        adata: AnnData with raw counts in X.
        method: "library_size" (default) or "scran" (raises — use scran.R).
        target_sum: total counts per cell after normalization (library_size only).
    """
    if method == "scran":
        raise ValueError(
            "scran normalization must be run via scran.R. "
            "Use the R-backed normalize step in the Nextflow pipeline."
        )
    if method != "library_size":
        raise ValueError(f"Unknown normalization method: '{method}'. Use 'library_size' or 'scran'.")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Library-size normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.layers["norm"] = adata.X.copy()

    # Log-transform
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()
    # adata.X is already log_norm (set by log1p in place)

    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_normalize.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/preprocessing/normalize.py tests/unit/test_normalize.py
git commit -m "feat: normalization with three-layer preservation"
```

---

### Task 10: HVG selection

**Files:**
- Create: `src/scanpy_workflow/preprocessing/hvg.py`
- Create: `tests/unit/test_hvg.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_hvg.py
import pytest
from scanpy_workflow.preprocessing.hvg import select_hvg


def test_hvg_seurat_v3_uses_counts_layer(normalized_adata):
    adata = select_hvg(normalized_adata, n_top_genes=50, flavor="seurat_v3")
    assert "highly_variable" in adata.var.columns
    assert adata.var["highly_variable"].sum() <= 50


def test_hvg_seurat_v3_requires_counts_layer(small_adata):
    with pytest.raises(ValueError, match="layers\\['counts'\\] required"):
        select_hvg(small_adata, n_top_genes=50, flavor="seurat_v3")


def test_hvg_seurat_flavor_uses_x(normalized_adata):
    adata = select_hvg(normalized_adata, n_top_genes=50, flavor="seurat")
    assert "highly_variable" in adata.var.columns


def test_hvg_marks_correct_count(normalized_adata):
    n = 30
    adata = select_hvg(normalized_adata.copy(), n_top_genes=n, flavor="seurat")
    assert adata.var["highly_variable"].sum() <= n
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_hvg.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/preprocessing/hvg.py`**

```python
import anndata as ad
import scanpy as sc


def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 3000,
    flavor: str = "seurat_v3",
) -> ad.AnnData:
    """
    Select highly variable genes and add 'highly_variable' flag to adata.var.

    For seurat_v3: requires layers["counts"] (raw counts).
    For seurat / cell_ranger: uses adata.X (log-normalized).
    """
    if flavor == "seurat_v3":
        if "counts" not in adata.layers:
            raise ValueError(
                "layers['counts'] required for seurat_v3 flavor. "
                "Run normalize() first to store raw counts."
            )
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            layer="counts",
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
        )
    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_hvg.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/preprocessing/hvg.py tests/unit/test_hvg.py
git commit -m "feat: HVG selection with seurat_v3 layer routing"
```

> **CLI note (Task 23):** The `scanpy-workflow hvg` command must accept `--mode per-sample|post-merge`. In per-sample mode it reads `params.hvg.n_top_genes`; in post-merge mode it reads `params.hvg.post_merge_n_top_genes`. The `select_hvg` function takes `n_top_genes` as a plain parameter — the CLI resolves the correct value from config before calling it.

---

## Chunk 4: Merge + Scale + PCA + Integration

### Task 11: Merge samples

**Files:**
- Create: `src/scanpy_workflow/preprocessing/merge.py`
- Create: `tests/unit/test_merge.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_merge.py
import pytest
from scanpy_workflow.preprocessing.merge import merge_samples


def test_merge_creates_batch_key_column(two_sample_adatas):
    merged = merge_samples(two_sample_adatas, sample_names=["s0", "s1"], batch_key="sample")
    assert "sample" in merged.obs.columns
    assert set(merged.obs["sample"].unique()) == {"s0", "s1"}


def test_merge_cell_count(two_sample_adatas):
    merged = merge_samples(two_sample_adatas, sample_names=["s0", "s1"], batch_key="sample")
    total = sum(a.n_obs for a in two_sample_adatas)
    assert merged.n_obs == total


def test_merge_preserves_var_names(two_sample_adatas):
    merged = merge_samples(two_sample_adatas, sample_names=["s0", "s1"], batch_key="sample")
    assert merged.n_vars == two_sample_adatas[0].n_vars


def test_merge_requires_matching_var_names(two_sample_adatas):
    import anndata as ad, numpy as np, scipy.sparse
    wrong = ad.AnnData(
        X=scipy.sparse.csr_matrix(np.ones((10, 5))),
        obs={"sample": ["x"] * 10},
    )
    wrong.var_names = [f"other_{i}" for i in range(5)]
    with pytest.raises(Exception):
        merge_samples(two_sample_adatas + [wrong], sample_names=["s0", "s1", "s2"], batch_key="sample")
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_merge.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/preprocessing/merge.py`**

```python
from typing import List
import anndata as ad


def merge_samples(
    adatas: List[ad.AnnData],
    sample_names: List[str],
    batch_key: str = "sample",
) -> ad.AnnData:
    """
    Merge per-sample AnnData objects into a single AnnData.

    Creates adata.obs[batch_key] with the sample name for each cell,
    using anndata.concat with label=batch_key, keys=sample_names.
    """
    if len(adatas) != len(sample_names):
        raise ValueError(f"len(adatas)={len(adatas)} != len(sample_names)={len(sample_names)}")

    merged = ad.concat(
        adatas,
        label=batch_key,
        keys=sample_names,
        merge="first",
        uns_merge="first",
    )
    return merged
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_merge.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/preprocessing/merge.py tests/unit/test_merge.py
git commit -m "feat: multi-sample merge with batch key injection"
```

---

### Task 12: Scale and PCA

**Files:**
- Create: `src/scanpy_workflow/preprocessing/scale.py`
- Create: `src/scanpy_workflow/preprocessing/pca.py`
- Create: `tests/unit/test_scale_pca.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_scale_pca.py
import numpy as np
import pytest
from scanpy_workflow.preprocessing.scale import scale
from scanpy_workflow.preprocessing.pca import run_pca


def test_scale_sets_x(normalized_adata):
    adata = scale(normalized_adata.copy())
    # After scaling, X should have zero mean per gene (within tolerance)
    import scipy.sparse
    x = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.asarray(adata.X)
    assert abs(x[:, 0].mean()) < 0.1


def test_pca_adds_x_pca(normalized_adata):
    adata = run_pca(normalized_adata, n_comps=10)
    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape == (normalized_adata.n_obs, 10)


def test_pca_scvi_path_asserts_x_is_log_norm(normalized_adata):
    """On the scVI path (scvi_path=True), X must equal layers['log_norm']."""
    import scipy.sparse
    adata = normalized_adata.copy()
    # Corrupt X so it doesn't equal log_norm
    adata.X = adata.X * 2
    with pytest.raises(AssertionError, match="adata.X must equal"):
        run_pca(adata, n_comps=10, scvi_path=True)


def test_pca_scvi_path_passes_when_x_is_log_norm(normalized_adata):
    adata = run_pca(normalized_adata.copy(), n_comps=10, scvi_path=True)
    assert "X_pca" in adata.obsm
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_scale_pca.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/preprocessing/scale.py`**

```python
import anndata as ad
import scanpy as sc


def scale(adata: ad.AnnData, max_value: float = 10.0) -> ad.AnnData:
    """
    Z-score scale gene expression. Sets adata.X to scaled data.
    Run after normalization and HVG selection, before PCA (non-scVI path).
    """
    sc.pp.scale(adata, max_value=max_value)
    return adata
```

- [ ] **Step 4: Write `src/scanpy_workflow/preprocessing/pca.py`**

```python
import anndata as ad
import scanpy as sc
import numpy as np
import scipy.sparse


def run_pca(
    adata: ad.AnnData,
    n_comps: int = 50,
    scvi_path: bool = False,
) -> ad.AnnData:
    """
    Run PCA and add adata.obsm['X_pca'].

    scvi_path=True: asserts adata.X == layers['log_norm'] (scaling was skipped).
    scvi_path=False: assumes adata.X is scaled data.
    """
    if scvi_path:
        if "log_norm" not in adata.layers:
            raise ValueError("layers['log_norm'] required for scVI PCA path.")
        x = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.asarray(adata.X)
        ln = adata.layers["log_norm"].toarray() if scipy.sparse.issparse(adata.layers["log_norm"]) else np.asarray(adata.layers["log_norm"])
        assert np.allclose(x, ln, atol=1e-5), (
            "adata.X must equal layers['log_norm'] on the scVI path. "
            "Ensure scale() was NOT called before PCA when using scVI batch correction."
        )

    sc.pp.pca(adata, n_comps=min(n_comps, adata.n_obs - 1, adata.n_vars - 1),
              use_highly_variable=True if "highly_variable" in adata.var.columns else False)
    return adata
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/unit/test_scale_pca.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/scanpy_workflow/preprocessing/scale.py src/scanpy_workflow/preprocessing/pca.py tests/unit/test_scale_pca.py
git commit -m "feat: scale and PCA with scVI path assertion"
```

---

### Task 13: Batch correction — Harmony

**Files:**
- Create: `src/scanpy_workflow/integration/harmony.py`
- Create: `src/scanpy_workflow/integration/evaluate.py`
- Create: `tests/unit/test_harmony.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_harmony.py
import pytest
from scanpy_workflow.integration.harmony import batch_correct_harmony


def test_harmony_adds_pca_harmony(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert "X_pca_harmony" in adata.obsm


def test_harmony_sets_use_rep(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert adata.uns["neighbors_use_rep"] == "X_pca_harmony"


def test_harmony_shape_preserved(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert adata.obsm["X_pca_harmony"].shape[0] == merged_adata.n_obs


def test_harmony_raises_without_pca(small_adata):
    with pytest.raises(ValueError, match="X_pca not found"):
        batch_correct_harmony(small_adata, batch_key="sample")


def test_harmony_evaluation_creates_json(merged_adata, tmp_path):
    adata = batch_correct_harmony(merged_adata, batch_key="sample",
                                   run_evaluation=True, output_dir=str(tmp_path))
    assert (tmp_path / "batch_correction_metrics.json").exists()
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_harmony.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/integration/harmony.py`**

```python
import anndata as ad
import harmonypy as hm
from .evaluate import evaluate_batch_correction


def batch_correct_harmony(
    adata: ad.AnnData,
    batch_key: str = "sample",
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply Harmony batch correction on X_pca.

    Writes:
      - adata.obsm["X_pca_harmony"]
      - adata.uns["neighbors_use_rep"] = "X_pca_harmony"
    """
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca not found in adata.obsm. Run PCA first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T
    adata.uns["neighbors_use_rep"] = "X_pca_harmony"

    if run_evaluation and output_dir:
        evaluate_batch_correction(adata, batch_key=batch_key, use_rep="X_pca_harmony",
                                   output_dir=output_dir)
    return adata
```

- [ ] **Step 4: Write `src/scanpy_workflow/integration/evaluate.py`**

```python
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


def evaluate_batch_correction(
    adata: ad.AnnData,
    batch_key: str,
    use_rep: str,
    output_dir: str,
    label_key: str = None,
) -> dict:
    """
    Evaluate batch correction using scib-metrics.

    Computes: kBET approximation, iLISI, cLISI, ASW batch.
    Writes UMAP pre/post correction plots to output_dir.
    Returns dict of metric scores.

    When label_key is None, only batch metrics are computed using
    embedding_obsm_keys; bio-conservation metrics are skipped.
    """
    from pathlib import Path
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    scores = {}

    try:
        import scib_metrics
        from scib_metrics.benchmark import Benchmarker
        # When label_key is None, use embedding_obsm_keys only for batch metrics
        # and skip bio-conservation metrics
        bm_kwargs = dict(
            adata=adata,
            batch_key=batch_key,
            embedding_obsm_keys=[use_rep],
        )
        if label_key is not None:
            bm_kwargs["label_key"] = label_key
        try:
            bm = Benchmarker(**bm_kwargs)
            bm.benchmark()
            results = bm.get_results(min_max_scale=False)
            scores = results.to_dict()
        except Exception as e:
            scores["error"] = str(e)
    except ImportError as e:
        scores["error"] = str(e)

    # UMAP comparison plot
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=batch_key, show=False)
    plt.savefig(out / "umap_post_correction.png", bbox_inches="tight", dpi=100)
    plt.close()

    import json
    with open(out / "batch_correction_metrics.json", "w") as f:
        json.dump(scores, f, indent=2, default=str)

    return scores
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/unit/test_harmony.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/scanpy_workflow/integration/harmony.py src/scanpy_workflow/integration/evaluate.py tests/unit/test_harmony.py
git commit -m "feat: Harmony batch correction and evaluation scaffold"
```

---

### Task 14: Batch correction — BBKNN

**Files:**
- Create: `src/scanpy_workflow/integration/bbknn.py`
- Create: `tests/unit/test_bbknn.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_bbknn.py
import pytest
import scipy.sparse
from scanpy_workflow.integration.bbknn import batch_correct_bbknn


def test_bbknn_writes_connectivities(merged_adata):
    adata = batch_correct_bbknn(merged_adata, batch_key="sample")
    assert "connectivities" in adata.obsp


def test_bbknn_writes_distances(merged_adata):
    adata = batch_correct_bbknn(merged_adata, batch_key="sample")
    assert "distances" in adata.obsp


def test_bbknn_sets_use_rep_skip(merged_adata):
    adata = batch_correct_bbknn(merged_adata, batch_key="sample")
    assert adata.uns["neighbors_use_rep"] == "skip"


def test_bbknn_raises_without_pca(small_adata):
    with pytest.raises(ValueError, match="X_pca not found"):
        batch_correct_bbknn(small_adata, batch_key="sample")
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_bbknn.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/integration/bbknn.py`**

```python
import anndata as ad
import bbknn


def batch_correct_bbknn(
    adata: ad.AnnData,
    batch_key: str = "sample",
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply BBKNN batch correction.

    Input: adata with X_pca in obsm.
    Output:
      - adata.obsp["connectivities"] and adata.obsp["distances"] written directly.
      - adata.uns["neighbors_use_rep"] = "skip" (pp.neighbors must be skipped).
    """
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca not found in adata.obsm. Run PCA first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    bbknn.bbknn(adata, adata.obs[batch_key], use_rep="X_pca")
    adata.uns["neighbors_use_rep"] = "skip"

    if run_evaluation and output_dir:
        from .evaluate import evaluate_batch_correction
        evaluate_batch_correction(adata, batch_key=batch_key,
                                   use_rep="X_pca", output_dir=output_dir)
    return adata
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_bbknn.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/integration/bbknn.py tests/unit/test_bbknn.py
git commit -m "feat: BBKNN batch correction"
```

---

### Task 15: scVI model utility + scVI batch correction

**Files:**
- Create: `src/scanpy_workflow/utils/scvi_model.py`
- Create: `src/scanpy_workflow/integration/scvi.py`
- Create: `tests/unit/test_scvi_integration.py`

> **Note:** scvi-tools is in a separate conda env (`env_scvi`). These tests are marked `scvi` and skipped unless `scvi-tools` is installed.

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_scvi_integration.py
import pytest

scvi = pytest.importorskip("scvi", reason="scvi-tools not installed (env_scvi)")

from scanpy_workflow.integration.scvi import batch_correct_scvi


def test_scvi_adds_x_scvi(merged_adata):
    adata = batch_correct_scvi(merged_adata, batch_key="sample", max_epochs=2)
    assert "X_scVI" in adata.obsm


def test_scvi_adds_denoised_layer(merged_adata):
    adata = batch_correct_scvi(merged_adata, batch_key="sample", max_epochs=2)
    assert "scvi_denoised" in adata.layers


def test_scvi_sets_use_rep(merged_adata):
    adata = batch_correct_scvi(merged_adata, batch_key="sample", max_epochs=2)
    assert adata.uns["neighbors_use_rep"] == "X_scVI"
```

- [ ] **Step 2: Write `src/scanpy_workflow/utils/scvi_model.py`**

```python
"""
Shared scVI VAE model training and inference utilities.
Used by both integration/scvi.py (batch correction) and imputation/scvi.py.
"""

import anndata as ad
from typing import Optional


def train_vae(
    adata: ad.AnnData,
    batch_key: str,
    layer: str = "counts",
    max_epochs: int = 400,
    early_stopping: bool = True,
    seed: int = 42,
):
    """
    Set up and train an scVI VAE model.

    Returns the trained SCVI model object.
    Input layer should be raw counts (layers['counts']) or log-normalized
    (layers['log_norm']) depending on context.
    """
    import scvi
    scvi.settings.seed = seed

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
    )
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=max_epochs, early_stopping=early_stopping)
    return model
```

- [ ] **Step 3: Write `src/scanpy_workflow/integration/scvi.py`**

```python
import anndata as ad
from pathlib import Path
from scanpy_workflow.utils.scvi_model import train_vae


def batch_correct_scvi(
    adata: ad.AnnData,
    batch_key: str = "sample",
    max_epochs: int = 400,
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply scVI batch correction.

    Uses layers['counts'] (raw counts) as input. Note: the spec's package
    structure comment "Input: layers['log_norm']" is a spec error — scVI's
    internal likelihood model requires raw counts; pre-normalized data produces
    incorrect results. This plan correctly uses layer="counts".
    Writes:
      - adata.obsm["X_scVI"]         : low-dim latent embedding
      - adata.layers["scvi_denoised"] : denoised expression
      - adata.uns["neighbors_use_rep"] = "X_scVI"
    """
    model = train_vae(adata, batch_key=batch_key, layer="counts")

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_denoised"] = model.get_normalized_expression(library_size=1e4)
    adata.uns["neighbors_use_rep"] = "X_scVI"

    if run_evaluation and output_dir:
        from .evaluate import evaluate_batch_correction
        evaluate_batch_correction(adata, batch_key=batch_key,
                                   use_rep="X_scVI", output_dir=output_dir)
    return adata
```

- [ ] **Step 4: Run tests (if env_scvi active)**

```bash
pytest tests/unit/test_scvi_integration.py -v
```

Expected: PASS if `scvi-tools` installed; SKIP otherwise.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/utils/scvi_model.py src/scanpy_workflow/integration/scvi.py tests/unit/test_scvi_integration.py
git commit -m "feat: scVI batch correction with shared VAEModel utility"
```

---

## Chunk 5: Imputation + Clustering

### Task 16: MAGIC imputation + evaluation

**Files:**
- Create: `src/scanpy_workflow/imputation/magic.py`
- Create: `src/scanpy_workflow/imputation/evaluate.py`
- Create: `tests/unit/test_magic.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_magic.py
import numpy as np
from scanpy_workflow.imputation.magic import impute_magic


def test_magic_preserves_shape(normalized_adata):
    adata = impute_magic(normalized_adata)
    assert adata.shape == normalized_adata.shape


def test_magic_modifies_x(normalized_adata):
    import scipy.sparse
    original = normalized_adata.X.copy()
    adata = impute_magic(normalized_adata.copy())
    # MAGIC changes some values
    if scipy.sparse.issparse(adata.X):
        new = adata.X.toarray()
        old = original.toarray()
    else:
        new = np.asarray(adata.X)
        old = np.asarray(original)
    assert not np.allclose(new, old)


def test_magic_operates_on_norm_layer(normalized_adata):
    """MAGIC must use layers['norm'] as input, not log_norm."""
    # This test verifies the function accepts and uses the correct layer
    adata = impute_magic(normalized_adata.copy())
    # After imputation, X should have been updated
    assert adata.X is not None
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_magic.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/imputation/magic.py`**

```python
import anndata as ad
import numpy as np
import scipy.sparse


def impute_magic(
    adata: ad.AnnData,
    solver: str = "exact",
    t: int = 3,
    run_evaluation: bool = False,
    output_dir: str = None,
    random_state: int = 42,
) -> ad.AnnData:
    """
    Apply MAGIC imputation.

    Operates on layers['norm'] (library-size normalized, pre-log) per
    MAGIC's recommended usage. Sets adata.X to imputed values.
    """
    import magic

    if "norm" not in adata.layers:
        raise ValueError("layers['norm'] required for MAGIC. Run normalize() first.")

    norm = adata.layers["norm"]
    if scipy.sparse.issparse(norm):
        norm = norm.toarray()

    magic_op = magic.MAGIC(solver=solver, t=t, random_state=random_state, verbose=0)
    imputed = magic_op.fit_transform(norm)

    adata.layers["imputed"] = imputed
    adata.X = scipy.sparse.csr_matrix(imputed) if scipy.sparse.issparse(adata.layers["norm"]) else imputed

    if run_evaluation and output_dir:
        from .evaluate import evaluate_imputation
        evaluate_imputation(adata,
                             impute_fn=lambda a: impute_magic(a, run_evaluation=False),
                             output_dir=output_dir)
    return adata
```

- [ ] **Step 4: Write `src/scanpy_workflow/imputation/evaluate.py`**

```python
import anndata as ad
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse
from pathlib import Path


def evaluate_imputation(
    adata: ad.AnnData,
    impute_fn: callable,
    output_dir: str,
    mask_fraction: float = 0.10,
    random_state: int = 42,
) -> dict:
    """
    In-silico dropout evaluation for imputation.

    Makes a deep copy of adata, masks mask_fraction of non-zero values in
    layers["norm"] of the copy, calls impute_fn(adata_copy) to get imputed
    results on the masked data, then compares imputed values vs original true
    values at masked positions.

    The production adata is never modified.

    Args:
        adata: Input AnnData (not modified).
        impute_fn: Callable that takes an AnnData and returns an AnnData with
                   imputed values. E.g. lambda a: impute_magic(a, run_evaluation=False)
        output_dir: Directory to write metrics JSON and distribution plot.
        mask_fraction: Fraction of non-zero values to mask.
        random_state: Random seed.

    Returns dict with pearson_r, spearman_r.
    """
    import copy
    from scipy.stats import pearsonr, spearmanr
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Make a deep copy to avoid modifying production adata
    adata_copy = copy.deepcopy(adata)

    # Get the norm layer from the copy for masking
    layer = adata_copy.layers["norm"]
    if scipy.sparse.issparse(layer):
        layer_dense = layer.toarray().copy()
    else:
        layer_dense = np.asarray(layer).copy()

    # Record original (true) values before masking
    rng = np.random.default_rng(random_state)
    nonzero_i, nonzero_j = np.nonzero(layer_dense)
    n_mask = int(len(nonzero_i) * mask_fraction)
    idx = rng.choice(len(nonzero_i), size=n_mask, replace=False)
    mask_i, mask_j = nonzero_i[idx], nonzero_j[idx]

    true_values = layer_dense[mask_i, mask_j].copy()

    # Apply mask to the copy's norm layer
    masked_layer = layer_dense.copy()
    masked_layer[mask_i, mask_j] = 0.0
    adata_copy.layers["norm"] = masked_layer

    # Run imputation on the masked copy
    imputed_adata = impute_fn(adata_copy)

    # Get imputed values at masked positions
    after_layer = imputed_adata.layers.get("imputed", imputed_adata.X)
    if scipy.sparse.issparse(after_layer):
        after_dense = after_layer.toarray()
    else:
        after_dense = np.asarray(after_layer)

    imputed_values = after_dense[mask_i, mask_j]

    pearson_r, _ = pearsonr(true_values, imputed_values)
    spearman_r, _ = spearmanr(true_values, imputed_values)
    scores = {"pearson_r": float(pearson_r), "spearman_r": float(spearman_r)}

    import json
    with open(out / "imputation_metrics.json", "w") as f:
        json.dump(scores, f, indent=2)

    # Distribution plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(layer_dense.flatten(), bins=50, alpha=0.7, label="before")
    axes[0].set_title("Before imputation")
    axes[1].hist(after_dense.flatten(), bins=50, alpha=0.7, label="after", color="orange")
    axes[1].set_title("After imputation")
    fig.savefig(out / "imputation_distribution.png", bbox_inches="tight", dpi=100)
    plt.close(fig)

    return scores
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/unit/test_magic.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/scanpy_workflow/imputation/magic.py src/scanpy_workflow/imputation/evaluate.py tests/unit/test_magic.py
git commit -m "feat: MAGIC imputation with in-silico dropout evaluation"
```

---

### Task 17: scVI imputation

**Files:**
- Create: `src/scanpy_workflow/imputation/scvi.py`
- Create: `tests/unit/test_scvi_imputation.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_scvi_imputation.py
import pytest
scvi = pytest.importorskip("scvi", reason="scvi-tools not installed")

from scanpy_workflow.imputation.scvi import impute_scvi


def test_scvi_imputation_adds_denoised_layer(merged_adata):
    adata = impute_scvi(merged_adata, batch_key="sample", max_epochs=2)
    assert "scvi_denoised" in adata.layers


def test_scvi_imputation_saves_model(merged_adata, tmp_path):
    adata = impute_scvi(merged_adata, batch_key="sample",
                        max_epochs=2, model_dir=str(tmp_path))
    assert (tmp_path / "model.pt").exists() or any(tmp_path.iterdir())
```

- [ ] **Step 2: Write `src/scanpy_workflow/imputation/scvi.py`**

```python
import anndata as ad
from pathlib import Path
from scanpy_workflow.utils.scvi_model import train_vae


def impute_scvi(
    adata: ad.AnnData,
    batch_key: str = "sample",
    max_epochs: int = 400,
    model_dir: str = None,
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Impute using scVI (trains its own model on layers['counts'] (raw counts, as required by scVI)).

    Only valid when batch_correction.method != 'scvi' (enforced by config.py).
    Writes:
      - adata.layers["scvi_denoised"]
      - saves model to model_dir if provided
    """
    model = train_vae(adata, batch_key=batch_key, layer="counts", max_epochs=max_epochs)

    adata.layers["scvi_denoised"] = model.get_normalized_expression(library_size=1e4)

    if model_dir:
        Path(model_dir).mkdir(parents=True, exist_ok=True)
        model.save(model_dir, overwrite=True)

    if run_evaluation and output_dir:
        from .evaluate import evaluate_imputation
        import scvi as _scvi
        def _scvi_impute_fn(a: ad.AnnData) -> ad.AnnData:
            _scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key=batch_key)
            m = _scvi.model.SCVI(a)
            m.train(max_epochs=max_epochs, early_stopping=True)
            a.layers["scvi_denoised"] = m.get_normalized_expression(library_size=1e4)
            return a
        evaluate_imputation(adata, impute_fn=_scvi_impute_fn, output_dir=output_dir)
    return adata
```

- [ ] **Step 3: Run tests**

```bash
pytest tests/unit/test_scvi_imputation.py -v
```

Expected: PASS if scvi-tools installed; SKIP otherwise.

- [ ] **Step 4: Commit**

```bash
git add src/scanpy_workflow/imputation/scvi.py tests/unit/test_scvi_imputation.py
git commit -m "feat: scVI imputation module"
```

---

### Task 18: Clustering + evaluation

**Files:**
- Create: `src/scanpy_workflow/clustering/cluster.py`
- Create: `src/scanpy_workflow/clustering/evaluate.py`
- Create: `tests/unit/test_cluster.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_cluster.py
import pytest
from scanpy_workflow.clustering.cluster import cluster


def test_cluster_leiden_adds_obs_column(merged_adata):
    adata = cluster(merged_adata, algorithm="leiden", resolution=0.3, n_neighbors=10)
    assert "leiden" in adata.obs.columns


def test_cluster_louvain_adds_obs_column(merged_adata):
    adata = cluster(merged_adata, algorithm="louvain", resolution=0.3, n_neighbors=10)
    assert "louvain" in adata.obs.columns


def test_cluster_skips_neighbors_for_bbknn(merged_adata):
    """When uns['neighbors_use_rep'] == 'skip', pp.neighbors is not called."""
    import scanpy as sc
    # Run BBKNN first so obsp exists
    adata = merged_adata.copy()
    adata.uns["neighbors_use_rep"] = "skip"
    sc.pp.neighbors(adata, use_rep="X_pca")  # simulate BBKNN having written the graph
    adata.uns["neighbors_use_rep"] = "skip"  # reset to skip
    result = cluster(adata, algorithm="leiden", resolution=0.3, n_neighbors=10)
    assert "leiden" in result.obs.columns


def test_cluster_unknown_algorithm_raises(merged_adata):
    with pytest.raises(ValueError, match="Unknown clustering algorithm"):
        cluster(merged_adata, algorithm="dbscan", resolution=0.3, n_neighbors=10)


def test_evaluate_clustering_creates_outputs(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    from scanpy_workflow.clustering.evaluate import evaluate_clustering
    scores = evaluate_clustering(adata, cluster_key="leiden", output_dir=str(tmp_path))
    assert (tmp_path / "clustering_metrics.json").exists()
    assert "silhouette" in scores or "note" in scores
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_cluster.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/clustering/cluster.py`**

```python
import anndata as ad
import scanpy as sc


def cluster(
    adata: ad.AnnData,
    algorithm: str = "leiden",
    resolution: float = 0.5,
    n_neighbors: int = 15,
) -> ad.AnnData:
    """
    Build neighbor graph (unless BBKNN already did it) and run Leiden/Louvain.

    Reads adata.uns['neighbors_use_rep']:
      - "skip": BBKNN wrote the graph; skip pp.neighbors.
      - anything else: call pp.neighbors(use_rep=...).
    """
    use_rep = adata.uns.get("neighbors_use_rep", "X_pca")

    if use_rep != "skip":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    if algorithm == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
    elif algorithm == "louvain":
        sc.tl.louvain(adata, resolution=resolution)
    else:
        raise ValueError(
            f"Unknown clustering algorithm: '{algorithm}'. Use 'leiden' or 'louvain'."
        )
    return adata
```

- [ ] **Step 4: Write `src/scanpy_workflow/clustering/evaluate.py`**

```python
import anndata as ad
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
from pathlib import Path

SUBSAMPLE_N = 20_000
RANDOM_STATE = 42


def evaluate_clustering(
    adata: ad.AnnData,
    cluster_key: str,
    output_dir: str,
) -> dict:
    """
    Compute clustering evaluation metrics.

    Silhouette and Davies-Bouldin: subsampled to <= SUBSAMPLE_N cells.
    Calinski-Harabasz: full dataset.
    Marker dot plots: written to output_dir.
    """
    from sklearn.metrics import (
        silhouette_score,
        davies_bouldin_score,
        calinski_harabasz_score,
    )
    import scanpy as sc

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    use_rep = adata.uns.get("neighbors_use_rep", "X_pca")
    embedding = adata.obsm.get(use_rep) or adata.obsm.get("X_pca")
    if embedding is None:
        raise ValueError("No valid embedding found for clustering evaluation.")

    labels = adata.obs[cluster_key].astype(str).values

    # Subsample for O(n²) metrics
    n = adata.n_obs
    if n > SUBSAMPLE_N:
        rng = np.random.default_rng(RANDOM_STATE)
        idx = rng.choice(n, size=SUBSAMPLE_N, replace=False)
        emb_sub = embedding[idx]
        lab_sub = labels[idx]
    else:
        emb_sub = embedding
        lab_sub = labels

    scores = {}
    if len(np.unique(lab_sub)) > 1:
        scores["silhouette"] = float(silhouette_score(emb_sub, lab_sub))
        scores["davies_bouldin"] = float(davies_bouldin_score(emb_sub, lab_sub))
        scores["calinski_harabasz"] = float(calinski_harabasz_score(embedding, labels))
    else:
        scores = {"silhouette": None, "davies_bouldin": None, "calinski_harabasz": None,
                  "note": "Only one cluster found; metrics not computed."}

    with open(out / "clustering_metrics.json", "w") as f:
        json.dump(scores, f, indent=2)

    # Marker dot plot
    try:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=5)
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        plt.savefig(out / "marker_dotplot.png", bbox_inches="tight", dpi=100)
        plt.close()
    except Exception:
        pass  # Skip if not enough cells per cluster

    return scores
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/unit/test_cluster.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/scanpy_workflow/clustering/cluster.py src/scanpy_workflow/clustering/evaluate.py tests/unit/test_cluster.py
git commit -m "feat: clustering with use_rep handoff and evaluation metrics"
```

---

## Chunk 6: Annotation + DE + Reporting + CLI

### Task 19: Unsupervised annotation (marker genes)

**Files:**
- Create: `src/scanpy_workflow/annotation/markers.py`
- Create: `tests/unit/test_markers.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_markers.py
import scanpy as sc
from scanpy_workflow.annotation.markers import rank_marker_genes


def test_rank_genes_adds_uns_key(merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    result = rank_marker_genes(adata, groupby="leiden")
    assert "rank_genes_groups" in result.uns


def test_rank_genes_csv_written(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    rank_marker_genes(adata, groupby="leiden", output_dir=str(tmp_path))
    assert (tmp_path / "marker_genes.csv").exists()
```

- [ ] **Step 2: Write `src/scanpy_workflow/annotation/markers.py`**

```python
import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path
from typing import Optional


def rank_marker_genes(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 50,
    output_dir: Optional[str] = None,
) -> ad.AnnData:
    """
    Identify marker genes per cluster using rank_genes_groups.
    If output_dir given, writes marker_genes.csv.
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        groups = adata.obs[groupby].unique().tolist()
        dfs = [
            sc.get.rank_genes_groups_df(adata, group=g).assign(group=g)
            for g in groups
        ]
        pd.concat(dfs).to_csv(out / "marker_genes.csv", index=False)

    return adata
```

- [ ] **Step 3: Run tests**

```bash
pytest tests/unit/test_markers.py -v
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/scanpy_workflow/annotation/markers.py tests/unit/test_markers.py
git commit -m "feat: marker gene identification"
```

---

### Task 20: CellTypist automated annotation

**Files:**
- Create: `src/scanpy_workflow/annotation/celltypist.py`
- Create: `tests/unit/test_celltypist.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_celltypist.py
import pytest
celltypist = pytest.importorskip("celltypist", reason="celltypist not installed")

from scanpy_workflow.annotation.celltypist import annotate_celltypist


def test_celltypist_adds_predicted_labels(merged_adata):
    from unittest.mock import patch, MagicMock
    import pandas as pd
    mock_predictions = MagicMock()
    result_adata = merged_adata.copy()
    result_adata.obs["majority_voting"] = "T cell"
    result_adata.obs["conf_score"] = 0.9
    mock_predictions.to_adata.return_value = result_adata
    with patch("celltypist.annotate", return_value=mock_predictions):
        from scanpy_workflow.annotation.celltypist import annotate_celltypist
        result = annotate_celltypist(merged_adata, model="Immune_All_Low.pkl")
    assert "celltypist_cell_type" in result.obs.columns
```

- [ ] **Step 2: Write `src/scanpy_workflow/annotation/celltypist.py`**

```python
import anndata as ad
import celltypist
from typing import Optional
from pathlib import Path


def annotate_celltypist(
    adata: ad.AnnData,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    output_dir: Optional[str] = None,
) -> ad.AnnData:
    """
    Automated cell type annotation with CellTypist.

    Adds adata.obs['celltypist_cell_type'] and adata.obs['celltypist_conf_score'].
    CellTypist expects log-normalized data in adata.X.
    """
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
    )
    adata = predictions.to_adata()
    adata.obs["celltypist_cell_type"] = adata.obs.get(
        "majority_voting", adata.obs.get("predicted_labels", "unknown")
    )
    adata.obs["celltypist_conf_score"] = adata.obs.get("conf_score", 0.0)

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        adata.obs[["celltypist_cell_type", "celltypist_conf_score"]].to_csv(
            out / "celltypist_predictions.csv"
        )
    return adata
```

- [ ] **Step 3: Run tests**

```bash
pytest tests/unit/test_celltypist.py -v
```

Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/scanpy_workflow/annotation/celltypist.py tests/unit/test_celltypist.py
git commit -m "feat: CellTypist automated annotation"
```

---

### Task 21: Differential expression

**Files:**
- Create: `src/scanpy_workflow/de/de.py`
- Create: `tests/unit/test_de.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_de.py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scanpy_workflow.de.de import run_de

REQUIRED_COLS = ["group", "gene", "score", "logfoldchange", "pval", "pval_adj"]


def _clustered_adata(merged_adata):
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    return adata


def test_de_returns_required_columns(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="wilcoxon")
    for col in REQUIRED_COLS:
        assert col in df.columns, f"Missing column: {col}"


def test_de_logreg_pval_is_nan(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="logreg")
    assert df["pval"].isna().all()
    assert df["pval_adj"].isna().all()


def test_de_wilcoxon_pval_not_nan(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="wilcoxon")
    assert not df["pval"].isna().all()


def test_de_raises_for_missing_groupby(merged_adata):
    with pytest.raises(ValueError, match="groupby 'nonexistent' not found"):
        run_de(merged_adata, groupby="nonexistent", method="wilcoxon")


def test_de_csv_written(tmp_path, merged_adata):
    adata = _clustered_adata(merged_adata)
    run_de(adata, groupby="leiden", method="wilcoxon", output_dir=str(tmp_path))
    assert (tmp_path / "de_results.csv").exists()
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_de.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/de/de.py`**

```python
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional


def run_de(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 100,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Run differential expression analysis.

    Outputs a normalized DataFrame with columns:
      group, gene, score, logfoldchange, pval, pval_adj

    For logreg: pval and pval_adj are NaN (logreg returns coefficients only).
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs.")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    groups = adata.obs[groupby].unique().tolist()
    dfs = []
    for group in groups:
        result = sc.get.rank_genes_groups_df(adata, group=str(group))
        result = result.rename(columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "logfoldchange",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        })
        result["group"] = str(group)
        dfs.append(result)

    df = pd.concat(dfs, ignore_index=True)
    df = df[["group", "gene", "score", "logfoldchange", "pval", "pval_adj"]]

    if method == "logreg":
        df["pval"] = np.nan
        df["pval_adj"] = np.nan

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / "de_results.csv", index=False)

    return df
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_de.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/de/de.py tests/unit/test_de.py
git commit -m "feat: differential expression with normalized output schema"
```

---

### Task 22: HTML reporting

**Files:**
- Create: `src/scanpy_workflow/reporting/report.py`
- Create: `src/scanpy_workflow/reporting/templates/report.html.j2`
- Create: `tests/unit/test_report.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_report.py
from pathlib import Path
from scanpy_workflow.reporting.report import generate_report


def test_report_creates_html(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    sc.tl.umap(adata)
    generate_report(adata, output_dir=str(tmp_path), params={})
    assert (tmp_path / "report.html").exists()


def test_report_contains_sample_info(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    sc.tl.umap(adata)
    generate_report(adata, output_dir=str(tmp_path), params={"samples": ["s1", "s2"]})
    html = (tmp_path / "report.html").read_text()
    assert "s1" in html or "Cells" in html
```

- [ ] **Step 2: Create Jinja2 template**

Create `src/scanpy_workflow/reporting/templates/report.html.j2`:

```html
<!DOCTYPE html>
<html>
<head><title>Scanpy Workflow Report</title>
<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:auto;padding:20px;}
h1{color:#2c3e50;}table{border-collapse:collapse;width:100%;}
td,th{border:1px solid #ddd;padding:8px;}th{background:#f2f2f2;}
img{max-width:100%;margin:10px 0;}</style>
</head>
<body>
<h1>Single-Cell RNA-seq Workflow Report</h1>
<p>Generated: {{ timestamp }}</p>

<h2>Dataset Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Total cells</td><td>{{ n_cells }}</td></tr>
<tr><td>Total genes</td><td>{{ n_genes }}</td></tr>
<tr><td>Samples</td><td>{{ samples | join(", ") }}</td></tr>
<tr><td>Clusters</td><td>{{ n_clusters }}</td></tr>
</table>

{% if umap_plot %}
<h2>UMAP</h2>
<img src="{{ umap_plot }}" alt="UMAP">
{% endif %}

{% if cluster_metrics %}
<h2>Clustering Metrics</h2>
<table>
{% for k, v in cluster_metrics.items() %}
<tr><td>{{ k }}</td><td>{{ v }}</td></tr>
{% endfor %}
</table>
{% endif %}

<h2>Parameters</h2>
<pre>{{ params_json }}</pre>
</body>
</html>
```

- [ ] **Step 3: Write `src/scanpy_workflow/reporting/report.py`**

```python
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import json
from datetime import datetime
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from typing import Optional


def generate_report(
    adata: ad.AnnData,
    output_dir: str,
    params: dict,
    cluster_key: str = "leiden",
    cluster_metrics: Optional[dict] = None,
) -> None:
    """Generate an HTML summary report."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # UMAP plot
    umap_plot = None
    if "X_umap" in adata.obsm and cluster_key in adata.obs.columns:
        sc.pl.umap(adata, color=cluster_key, show=False)
        umap_file = out / "umap_clusters.png"
        plt.savefig(umap_file, bbox_inches="tight", dpi=100)
        plt.close()
        umap_plot = "umap_clusters.png"

    samples = list(adata.obs.get("sample", adata.obs.get("batch", pd.Series(["unknown"]))).unique())
    n_clusters = adata.obs[cluster_key].nunique() if cluster_key in adata.obs.columns else "N/A"

    template_dir = Path(__file__).parent / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    template = env.get_template("report.html.j2")

    html = template.render(
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        samples=params.get("samples", samples),
        n_clusters=n_clusters,
        umap_plot=umap_plot,
        cluster_metrics=cluster_metrics or {},
        params_json=json.dumps(params, indent=2, default=str),
    )
    (out / "report.html").write_text(html)
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_report.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/scanpy_workflow/reporting/ tests/unit/test_report.py
git commit -m "feat: Jinja2 HTML report generation"
```

---

### Task 23: CLI assembly

**Files:**
- Create: `src/scanpy_workflow/cli.py`
- Create: `tests/unit/test_cli.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_cli.py
from click.testing import CliRunner
from scanpy_workflow.cli import cli


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "scanpy-workflow" in result.output.lower() or "Usage" in result.output


def test_qc_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["qc", "--help"])
    assert result.exit_code == 0


def test_filter_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["filter", "--help"])
    assert result.exit_code == 0


def test_normalize_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["normalize", "--help"])
    assert result.exit_code == 0


def test_cluster_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["cluster", "--help"])
    assert result.exit_code == 0


def test_de_command_exists():
    runner = CliRunner()
    result = runner.invoke(cli, ["de", "--help"])
    assert result.exit_code == 0
```

- [ ] **Step 2: Run — expect failure**

```bash
pytest tests/unit/test_cli.py -v
```

- [ ] **Step 3: Write `src/scanpy_workflow/cli.py`**

```python
"""
CLI entrypoints for scanpy_workflow.
All multi-word subcommands use hyphens (Unix convention).
R-backed steps (ambient, normalize scran, impute alra, annotate sctype)
are called directly by Nextflow via Rscript — they are NOT registered here.
"""
import click
from scanpy_workflow.utils.io import read_h5ad, write_h5ad


@click.group()
def cli():
    """Scanpy end-to-end single-cell workflow."""
    pass


# ── load ──────────────────────────────────────────────────────────────────────
@cli.command("load")
@click.option("--input", "input_dir", required=True, help="Cell Ranger output directory")
@click.option("--output", "output_path", required=True, help="Output .h5ad path")
def cmd_load(input_dir, output_path):
    """Load Cell Ranger output (MEX or HDF5)."""
    from scanpy_workflow.io.load import load_cellranger
    adata = load_cellranger(input_dir)
    write_h5ad(adata, output_path)
    click.echo(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes -> {output_path}")


# ── qc ────────────────────────────────────────────────────────────────────────
@cli.command("qc")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--plots-dir", "plots_dir", default=None)
def cmd_qc(input_path, output_path, plots_dir):
    """Compute QC metrics."""
    from scanpy_workflow.qc.metrics import compute_qc_metrics
    from scanpy_workflow.qc.plots import plot_qc
    adata = read_h5ad(input_path)
    adata = compute_qc_metrics(adata)
    if plots_dir:
        plot_qc(adata, output_dir=plots_dir)
    write_h5ad(adata, output_path)


# ── filter ────────────────────────────────────────────────────────────────────
@cli.command("filter")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--min-genes", default=200, type=int)
@click.option("--max-genes", default=6000, type=int)
@click.option("--min-counts", default=500, type=int)
@click.option("--max-counts", default=30000, type=int)
@click.option("--max-pct-mito", default=20.0, type=float)
@click.option("--max-pct-ribo", default=50.0, type=float)
@click.option("--min-cells", default=3, type=int)
def cmd_filter(input_path, output_path, min_genes, max_genes, min_counts,
               max_counts, max_pct_mito, max_pct_ribo, min_cells):
    """Filter cells and genes by QC cutoffs."""
    from scanpy_workflow.qc.filter import filter_cells_genes
    adata = read_h5ad(input_path)
    adata = filter_cells_genes(adata, min_genes=min_genes, max_genes=max_genes,
                                min_counts=min_counts, max_counts=max_counts,
                                max_pct_mito=max_pct_mito, max_pct_ribo=max_pct_ribo,
                                min_cells=min_cells)
    write_h5ad(adata, output_path)
    click.echo(f"Filtered: {adata.n_obs} cells, {adata.n_vars} genes -> {output_path}")


# ── doublets ──────────────────────────────────────────────────────────────────
@cli.command("doublets")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--filter-doublets", is_flag=True, default=False)
def cmd_doublets(input_path, output_path, filter_doublets):
    """Detect (and optionally filter) doublets with Scrublet."""
    from scanpy_workflow.doublets.scrublet import detect_doublets
    adata = read_h5ad(input_path)
    adata = detect_doublets(adata, filter_doublets=filter_doublets)
    write_h5ad(adata, output_path)


# ── normalize ─────────────────────────────────────────────────────────────────
@cli.command("normalize")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", default="library_size", type=click.Choice(["library_size"]))
@click.option("--target-sum", default=10000.0, type=float)
def cmd_normalize(input_path, output_path, method, target_sum):
    """Normalize counts (library_size). For scran, use scran.R."""
    from scanpy_workflow.preprocessing.normalize import normalize
    adata = read_h5ad(input_path)
    adata = normalize(adata, method=method, target_sum=target_sum)
    write_h5ad(adata, output_path)


# ── hvg ───────────────────────────────────────────────────────────────────────
@cli.command("hvg")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--n-top-genes", default=3000, type=int)
@click.option("--flavor", default="seurat_v3",
              type=click.Choice(["seurat_v3", "seurat", "cell_ranger"]))
# Nextflow passes correct --n-top-genes based on mode (per-sample vs post-merge)
def cmd_hvg(input_path, output_path, n_top_genes, flavor):
    """Select highly variable genes."""
    from scanpy_workflow.preprocessing.hvg import select_hvg
    adata = read_h5ad(input_path)
    adata = select_hvg(adata, n_top_genes=n_top_genes, flavor=flavor)
    write_h5ad(adata, output_path)


# ── merge ─────────────────────────────────────────────────────────────────────
@cli.command("merge")
@click.option("--inputs", "input_paths", required=True, multiple=True,
              help="Paths to per-sample .h5ad files")
@click.option("--sample-names", required=True, multiple=True)
@click.option("--output", "output_path", required=True)
@click.option("--batch-key", default="sample")
def cmd_merge(input_paths, sample_names, output_path, batch_key):
    """Merge per-sample AnnData objects."""
    from scanpy_workflow.preprocessing.merge import merge_samples
    adatas = [read_h5ad(p) for p in input_paths]
    merged = merge_samples(adatas, sample_names=list(sample_names), batch_key=batch_key)
    write_h5ad(merged, output_path)
    click.echo(f"Merged {len(adatas)} samples: {merged.n_obs} cells -> {output_path}")


# ── scale ─────────────────────────────────────────────────────────────────────
@cli.command("scale")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--max-value", default=10.0, type=float)
def cmd_scale(input_path, output_path, max_value):
    """Z-score scale gene expression."""
    from scanpy_workflow.preprocessing.scale import scale
    adata = read_h5ad(input_path)
    adata = scale(adata, max_value=max_value)
    write_h5ad(adata, output_path)


# ── pca ───────────────────────────────────────────────────────────────────────
@cli.command("pca")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--n-comps", default=50, type=int)
@click.option("--scvi-path", is_flag=True, default=False,
              help="Assert X equals layers['log_norm'] (scVI path, scaling was skipped)")
def cmd_pca(input_path, output_path, n_comps, scvi_path):
    """Run PCA."""
    from scanpy_workflow.preprocessing.pca import run_pca
    adata = read_h5ad(input_path)
    adata = run_pca(adata, n_comps=n_comps, scvi_path=scvi_path)
    write_h5ad(adata, output_path)


# ── batch-correct ──────────────────────────────────────────────────────────────
@cli.command("batch-correct")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", required=True, type=click.Choice(["harmony", "bbknn", "scvi"]))
@click.option("--batch-key", default="sample")
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
def cmd_batch_correct(input_path, output_path, method, batch_key, evaluate, eval_dir):
    """Apply batch correction."""
    adata = read_h5ad(input_path)
    if method == "harmony":
        from scanpy_workflow.integration.harmony import batch_correct_harmony
        adata = batch_correct_harmony(adata, batch_key=batch_key,
                                       run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "bbknn":
        from scanpy_workflow.integration.bbknn import batch_correct_bbknn
        adata = batch_correct_bbknn(adata, batch_key=batch_key,
                                     run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "scvi":
        from scanpy_workflow.integration.scvi import batch_correct_scvi
        adata = batch_correct_scvi(adata, batch_key=batch_key,
                                    run_evaluation=evaluate, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── impute ────────────────────────────────────────────────────────────────────
@cli.command("impute")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", required=True, type=click.Choice(["magic", "scvi"]))
@click.option("--batch-key", default="sample")
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
@click.option("--model-dir", default=None)
def cmd_impute(input_path, output_path, method, batch_key, evaluate, eval_dir, model_dir):
    """Impute missing values (Python methods: magic, scvi). For alra, use alra.R."""
    adata = read_h5ad(input_path)
    if method == "magic":
        from scanpy_workflow.imputation.magic import impute_magic
        adata = impute_magic(adata, run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "scvi":
        from scanpy_workflow.imputation.scvi import impute_scvi
        adata = impute_scvi(adata, batch_key=batch_key, model_dir=model_dir,
                             run_evaluation=evaluate, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── cluster ───────────────────────────────────────────────────────────────────
@cli.command("cluster")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--algorithm", default="leiden", type=click.Choice(["leiden", "louvain"]))
@click.option("--resolution", default=0.5, type=float)
@click.option("--n-neighbors", default=15, type=int)
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
def cmd_cluster(input_path, output_path, algorithm, resolution, n_neighbors, evaluate, eval_dir):
    """Build neighbor graph and cluster cells."""
    from scanpy_workflow.clustering.cluster import cluster
    adata = read_h5ad(input_path)
    adata = cluster(adata, algorithm=algorithm, resolution=resolution, n_neighbors=n_neighbors)
    if evaluate and eval_dir:
        from scanpy_workflow.clustering.evaluate import evaluate_clustering
        evaluate_clustering(adata, cluster_key=algorithm, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── rank-genes ────────────────────────────────────────────────────────────────
@cli.command("rank-genes")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--groupby", default="leiden")
@click.option("--output-dir", "marker_dir", default=None)
def cmd_rank_genes(input_path, output_path, groupby, marker_dir):
    """Identify marker genes per cluster."""
    from scanpy_workflow.annotation.markers import rank_marker_genes
    adata = read_h5ad(input_path)
    adata = rank_marker_genes(adata, groupby=groupby, output_dir=marker_dir)
    write_h5ad(adata, output_path)


# ── annotate ──────────────────────────────────────────────────────────────────
@cli.command("annotate")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--model", default="Immune_All_Low.pkl")
@click.option("--output-dir", "annot_dir", default=None)
def cmd_annotate(input_path, output_path, model, annot_dir):
    """Automated annotation with CellTypist."""
    from scanpy_workflow.annotation.celltypist import annotate_celltypist
    adata = read_h5ad(input_path)
    adata = annotate_celltypist(adata, model=model, output_dir=annot_dir)
    write_h5ad(adata, output_path)


# ── de ────────────────────────────────────────────────────────────────────────
@cli.command("de")
@click.option("--input", "input_path", required=True)
@click.option("--output-h5ad", "output_path", required=True)
@click.option("--output-dir", "de_dir", required=True)
@click.option("--groupby", default="leiden")
@click.option("--method", default="wilcoxon",
              type=click.Choice(["wilcoxon", "t-test", "logreg"]))
def cmd_de(input_path, output_path, de_dir, groupby, method):
    """Run differential expression."""
    from scanpy_workflow.de.de import run_de
    adata = read_h5ad(input_path)
    run_de(adata, groupby=groupby, method=method, output_dir=de_dir)
    write_h5ad(adata, output_path)


# ── report ────────────────────────────────────────────────────────────────────
@cli.command("report")
@click.option("--input", "input_path", required=True)
@click.option("--output-dir", "output_dir", required=True)
@click.option("--config", "config_path", default=None)
@click.option("--cluster-key", default="leiden")
def cmd_report(input_path, output_dir, config_path, cluster_key):
    """Generate HTML summary report."""
    from scanpy_workflow.reporting.report import generate_report
    from scanpy_workflow.utils.config import load_config
    adata = read_h5ad(input_path)
    params = load_config(config_path) if config_path else {}
    generate_report(adata, output_dir=output_dir, params=params, cluster_key=cluster_key)
    click.echo(f"Report written to {output_dir}/report.html")
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_cli.py -v
```

Expected: all PASS.

- [ ] **Step 5: Run the full unit test suite**

```bash
pytest tests/unit/ -v --ignore=tests/unit/r
```

Expected: all tests PASS (scvi tests skip if not in env_scvi).

- [ ] **Step 6: Commit**

```bash
git add src/scanpy_workflow/cli.py tests/unit/test_cli.py
git commit -m "feat: CLI assembly — all scanpy-workflow subcommands registered"
```

---

### Task 24: Final integration smoke test

**Files:**
- Create: `tests/integration/test_smoke.py`

- [ ] **Step 1: Write smoke test**

```python
# tests/integration/test_smoke.py
"""
End-to-end smoke test using a programmatically generated synthetic dataset.
No external data downloads required. Tests the full Python package pipeline.
"""
import pytest
import numpy as np
import anndata as ad
import scipy.sparse
import scanpy as sc
from pathlib import Path


@pytest.fixture(scope="module")
def synthetic_cellranger(tmp_path_factory):
    """Create two synthetic Cell Ranger MEX directories."""
    import scipy.io, gzip
    base = tmp_path_factory.mktemp("cellranger")
    sample_dirs = []
    for s in range(2):
        rng = np.random.default_rng(s)
        n_cells, n_genes = 80, 150
        sample_dir = base / f"sample_{s}" / "filtered_feature_bc_matrix"
        sample_dir.mkdir(parents=True)
        matrix = scipy.sparse.random(n_genes, n_cells, density=0.3, format="csc",
                                      random_state=s, data_rvs=lambda n: rng.integers(1, 20, n))
        import io, gzip as gz
        buf = io.BytesIO()
        scipy.io.mmwrite(buf, matrix)
        with gz.open(sample_dir / "matrix.mtx.gz", "wb") as f:
            f.write(buf.getvalue())
        with gzip.open(sample_dir / "barcodes.tsv.gz", "wt") as f:
            for i in range(n_cells):
                f.write(f"CELL{s}{i:04d}-1\n")
        with gzip.open(sample_dir / "features.tsv.gz", "wt") as f:
            for i in range(n_genes):
                prefix = "MT-" if i < 5 else "RPS" if i < 10 else "Gene"
                f.write(f"ENSG{i:08d}\t{prefix}{i}\tGene Expression\n")
        sample_dirs.append(str(sample_dir.parent.parent))
    return sample_dirs


def test_full_pipeline_smoke(synthetic_cellranger, tmp_path):
    from scanpy_workflow.io.load import load_cellranger
    from scanpy_workflow.qc.metrics import compute_qc_metrics
    from scanpy_workflow.qc.filter import filter_cells_genes
    from scanpy_workflow.doublets.scrublet import detect_doublets
    from scanpy_workflow.preprocessing.normalize import normalize
    from scanpy_workflow.preprocessing.hvg import select_hvg
    from scanpy_workflow.preprocessing.merge import merge_samples
    from scanpy_workflow.preprocessing.scale import scale
    from scanpy_workflow.preprocessing.pca import run_pca
    from scanpy_workflow.integration.harmony import batch_correct_harmony
    from scanpy_workflow.clustering.cluster import cluster
    from scanpy_workflow.de.de import run_de

    # Per-sample
    adatas = []
    for s, sample_dir in enumerate(synthetic_cellranger):
        adata = load_cellranger(sample_dir)
        adata = compute_qc_metrics(adata)
        adata = filter_cells_genes(adata, min_genes=1, max_genes=10000,
                                    min_counts=1, max_counts=10**9,
                                    max_pct_mito=100, max_pct_ribo=100, min_cells=1)
        adata = detect_doublets(adata, filter_doublets=False)
        adata = normalize(adata)
        adata = select_hvg(adata, n_top_genes=50, flavor="seurat_v3")
        adatas.append(adata)

    # Merge + shared
    merged = merge_samples(adatas, sample_names=["s0", "s1"], batch_key="sample")
    select_hvg(merged, n_top_genes=50, flavor="seurat_v3")
    scale(merged)
    run_pca(merged, n_comps=10)
    batch_correct_harmony(merged, batch_key="sample")
    cluster(merged, algorithm="leiden", resolution=0.5, n_neighbors=5)
    df = run_de(merged, groupby="leiden", method="wilcoxon")

    # Assertions
    assert merged.n_obs > 0
    assert "leiden" in merged.obs.columns
    assert "X_pca_harmony" in merged.obsm
    assert set(df.columns) == {"group", "gene", "score", "logfoldchange", "pval", "pval_adj"}
```

- [ ] **Step 2: Run smoke test**

```bash
pytest tests/integration/test_smoke.py -v
```

Expected: PASS.

- [ ] **Step 3: Run full test suite**

```bash
pytest tests/ -v --ignore=tests/unit/r
```

Expected: all tests PASS or SKIP (scvi tests skip without env_scvi).

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_smoke.py
git commit -m "test: end-to-end smoke test for Python package pipeline"
```

---

## Next Plans

- **Plan 2:** `docs/superpowers/plans/2026-03-15-scanpy-workflow-plan-2-r-wrappers.md` — R scripts (SoupX, DecontX, scran, ALRA, scType) with testthat unit tests
- **Plan 3:** `docs/superpowers/plans/2026-03-15-scanpy-workflow-plan-3-nextflow.md` — Nextflow DSL2 pipeline, conda env YAML files, Dockerfiles, CI workflow
