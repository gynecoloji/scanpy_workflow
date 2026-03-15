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
