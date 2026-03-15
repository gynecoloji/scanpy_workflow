import numpy as np
from scanpy_workflow.doublets.scrublet import detect_doublets


def test_scrublet_adds_doublet_score(small_adata):
    adata = detect_doublets(small_adata)
    assert "doublet_score" in adata.obs.columns
    assert "predicted_doublet" in adata.obs.columns


def test_scrublet_doublet_score_range(small_adata):
    adata = detect_doublets(small_adata)
    assert (adata.obs["doublet_score"] >= 0).all()
    assert (adata.obs["doublet_score"] <= 1).all()


def test_scrublet_filter_removes_doublets(small_adata):
    adata = detect_doublets(small_adata, filter_doublets=True)
    # After filtering, no predicted doublets remain
    assert not adata.obs["predicted_doublet"].any()
