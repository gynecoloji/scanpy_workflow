import pytest
from scanpy_workflow.integration.harmony import batch_correct_harmony


def test_harmony_adds_pca_harmony(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert "X_pca_harmony" in adata.obsm


def test_harmony_sets_use_rep(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert adata.uns["neighbors_use_rep"] == "X_pca_harmony"


def test_harmony_shape_preserved(merged_adata):
    adata = batch_correct_harmony(merged_adata, batch_key="sample")
    assert adata.obsm["X_pca_harmony"].shape[0] == merged_adata.n_obs


def test_harmony_raises_without_pca(small_adata):
    with pytest.raises(ValueError, match="X_pca not found"):
        batch_correct_harmony(small_adata, batch_key="sample")


def test_harmony_evaluation_creates_json(merged_adata, tmp_path):
    adata = batch_correct_harmony(merged_adata, batch_key="sample",
                                   run_evaluation=True, output_dir=str(tmp_path))
    assert (tmp_path / "batch_correction_metrics.json").exists()
