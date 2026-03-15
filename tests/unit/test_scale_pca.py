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
