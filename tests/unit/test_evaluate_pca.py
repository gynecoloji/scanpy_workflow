import json
import numpy as np
import pytest
import scanpy as sc
import anndata as ad
import scipy.sparse

from scanpy_workflow.preprocessing.evaluate_pca import (
    evaluate_pca,
    _detect_elbow,
    _comps_for_threshold,
)


@pytest.fixture
def pca_adata(merged_adata):
    """merged_adata already has X_pca; we just need adata.uns['pca']."""
    adata = merged_adata.copy()
    # merged_adata has sc.pp.pca run with n_comps=20 in conftest
    return adata


# ── evaluate_pca ──────────────────────────────────────────────────────────────

def test_evaluate_pca_writes_json(tmp_path, pca_adata):
    scores = evaluate_pca(pca_adata, output_dir=str(tmp_path))
    assert (tmp_path / "pca_metrics.json").exists()
    with open(tmp_path / "pca_metrics.json") as f:
        data = json.load(f)
    assert "n_comps_used" in data
    assert "elbow_pc" in data
    assert "cumvar_80pct" in data
    assert "cumvar_90pct" in data
    assert "variance_explained" in data


def test_evaluate_pca_writes_scree_plot(tmp_path, pca_adata):
    evaluate_pca(pca_adata, output_dir=str(tmp_path))
    assert (tmp_path / "pca_scree.png").exists()


def test_evaluate_pca_returns_dict(tmp_path, pca_adata):
    scores = evaluate_pca(pca_adata, output_dir=str(tmp_path))
    assert isinstance(scores, dict)
    assert scores["n_comps_used"] > 0


def test_evaluate_pca_n_comps_matches_adata(tmp_path, pca_adata):
    scores = evaluate_pca(pca_adata, output_dir=str(tmp_path))
    n_pcs = len(pca_adata.uns["pca"]["variance_ratio"])
    assert scores["n_comps_used"] == n_pcs
    assert len(scores["variance_explained"]) == n_pcs


def test_evaluate_pca_cumvar_80_leq_90(tmp_path, pca_adata):
    scores = evaluate_pca(pca_adata, output_dir=str(tmp_path))
    assert scores["cumvar_80pct"] <= scores["cumvar_90pct"]


def test_evaluate_pca_elbow_in_range(tmp_path, pca_adata):
    scores = evaluate_pca(pca_adata, output_dir=str(tmp_path))
    n = scores["n_comps_used"]
    assert 1 <= scores["elbow_pc"] <= n


def test_evaluate_pca_raises_without_pca(tmp_path, small_adata):
    """Must raise if PCA has not been run."""
    with pytest.raises(ValueError, match="PCA results not found"):
        evaluate_pca(small_adata, output_dir=str(tmp_path))


def test_evaluate_pca_creates_output_dir(tmp_path, pca_adata):
    nested = tmp_path / "nested" / "eval"
    evaluate_pca(pca_adata, output_dir=str(nested))
    assert nested.exists()


# ── _detect_elbow ─────────────────────────────────────────────────────────────

def test_detect_elbow_returns_int():
    var = np.array([0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
    result = _detect_elbow(var)
    assert isinstance(result, int)
    assert 1 <= result <= len(var)


def test_detect_elbow_monotone_decreasing():
    """Elbow of a sharply decaying curve should be at the early PCs."""
    var = np.exp(-np.arange(10) * 0.8)
    elbow = _detect_elbow(var)
    assert elbow <= 4


def test_detect_elbow_short_array():
    assert _detect_elbow(np.array([0.5, 0.3])) in [1, 2]


# ── _comps_for_threshold ──────────────────────────────────────────────────────

def test_comps_for_80pct():
    cumvar = np.array([0.3, 0.55, 0.72, 0.84, 0.91, 0.96])
    n = _comps_for_threshold(cumvar, 0.80)
    assert n == 4   # index 3 (0.84) is first to cross 0.80


def test_comps_for_threshold_never_exceeded():
    """If threshold is higher than max cumvar, return the last PC."""
    cumvar = np.array([0.3, 0.55, 0.72])
    n = _comps_for_threshold(cumvar, 0.99)
    assert n == 3
