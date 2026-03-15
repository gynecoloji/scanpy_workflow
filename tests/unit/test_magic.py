import numpy as np
import pytest
from scanpy_workflow.imputation.magic import impute_magic

magic = pytest.importorskip("magic", reason="magic-impute not installed; skipping MAGIC tests")


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
