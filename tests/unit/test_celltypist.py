import pytest
celltypist = pytest.importorskip("celltypist", reason="celltypist not installed")

from scanpy_workflow.annotation.celltypist import annotate_celltypist


def test_celltypist_adds_predicted_labels(merged_adata):
    from unittest.mock import patch, MagicMock
    import pandas as pd
    mock_predictions = MagicMock()
    result_adata = merged_adata.copy()
    result_adata.obs["majority_voting"] = "T cell"
    result_adata.obs["conf_score"] = 0.9
    mock_predictions.to_adata.return_value = result_adata
    with patch("celltypist.annotate", return_value=mock_predictions):
        from scanpy_workflow.annotation.celltypist import annotate_celltypist
        result = annotate_celltypist(merged_adata, model="Immune_All_Low.pkl")
    assert "celltypist_cell_type" in result.obs.columns
