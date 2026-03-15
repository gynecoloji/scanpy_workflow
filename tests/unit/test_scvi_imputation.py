import pytest
scvi = pytest.importorskip("scvi", reason="scvi-tools not installed")

from scanpy_workflow.imputation.scvi import impute_scvi


def test_scvi_imputation_adds_denoised_layer(merged_adata):
    adata = impute_scvi(merged_adata, batch_key="sample", max_epochs=2)
    assert "scvi_denoised" in adata.layers


def test_scvi_imputation_saves_model(merged_adata, tmp_path):
    adata = impute_scvi(merged_adata, batch_key="sample",
                        max_epochs=2, model_dir=str(tmp_path))
    assert (tmp_path / "model.pt").exists() or any(tmp_path.iterdir())
