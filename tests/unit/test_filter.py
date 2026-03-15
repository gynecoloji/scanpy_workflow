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
