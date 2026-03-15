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
