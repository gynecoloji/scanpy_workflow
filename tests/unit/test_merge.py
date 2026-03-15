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
