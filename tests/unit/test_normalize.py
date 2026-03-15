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
