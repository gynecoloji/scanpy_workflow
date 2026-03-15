import pytest
import sys
import numpy as np
import anndata as ad
import scanpy as sc
import scipy.sparse
import gzip
import io
import os
import scipy.io
from pathlib import Path

# Add src directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


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
    adata = anndata.concat(
        two_sample_adatas,
        label="sample",
        keys=["sample_0", "sample_1"],
        merge="unique",
    )
    sc.pp.highly_variable_genes(adata, n_top_genes=100, flavor="seurat")
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=20, use_highly_variable=True)
    adata.uns["neighbors_use_rep"] = "X_pca"
    return adata


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
                         scipy.sparse.csc_matrix(rng.poisson(0.1, size=(n_genes, len(bc_list))).astype("float32")))
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
