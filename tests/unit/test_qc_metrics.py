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
