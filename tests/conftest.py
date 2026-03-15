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
