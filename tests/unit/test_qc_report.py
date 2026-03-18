import base64
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

from scanpy_workflow.qc.metrics import compute_qc_metrics
from scanpy_workflow.qc.report import generate_qc_report, _compute_filter_steps, _compute_metric_stats
from scanpy_workflow.qc.plots import (
    plot_violin_with_thresholds,
    plot_scatter_with_thresholds,
    plot_qc_histograms,
    plot_knee_curve,
    plot_filter_waterfall,
)

THRESHOLDS = dict(
    min_genes=100, max_genes=500,
    min_counts=200, max_counts=5000,
    max_pct_mito=20.0, max_pct_ribo=50.0,
)


@pytest.fixture
def qc_adata(small_adata):
    return compute_qc_metrics(small_adata.copy())


# ── report generation ─────────────────────────────────────────────────────────

def test_generate_qc_report_creates_html(tmp_path, qc_adata):
    out = tmp_path / "qc_report.html"
    generate_qc_report(qc_adata, sample_name="testSample",
                       thresholds=THRESHOLDS, output_path=str(out))
    assert out.exists()
    assert out.stat().st_size > 5000


def test_report_contains_sample_name(tmp_path, qc_adata):
    out = tmp_path / "qc_report.html"
    generate_qc_report(qc_adata, sample_name="MySample42",
                       thresholds=THRESHOLDS, output_path=str(out))
    html = out.read_text()
    assert "MySample42" in html


def test_report_contains_cell_counts(tmp_path, qc_adata):
    out = tmp_path / "qc_report.html"
    generate_qc_report(qc_adata, sample_name="s",
                       thresholds=THRESHOLDS, output_path=str(out))
    html = out.read_text()
    # The total cell count should appear somewhere in the report
    assert str(qc_adata.n_obs) in html


def test_report_contains_embedded_images(tmp_path, qc_adata):
    out = tmp_path / "qc_report.html"
    generate_qc_report(qc_adata, sample_name="s",
                       thresholds=THRESHOLDS, output_path=str(out))
    html = out.read_text()
    assert "data:image/png;base64," in html


def test_report_raises_without_qc_metrics(tmp_path, small_adata):
    """Report must raise if compute_qc_metrics has not been run."""
    out = tmp_path / "qc_report.html"
    with pytest.raises(ValueError, match="Missing QC columns"):
        generate_qc_report(small_adata, sample_name="s",
                           thresholds=THRESHOLDS, output_path=str(out))


def test_report_creates_parent_dirs(tmp_path, qc_adata):
    out = tmp_path / "nested" / "deep" / "qc_report.html"
    generate_qc_report(qc_adata, sample_name="s",
                       thresholds=THRESHOLDS, output_path=str(out))
    assert out.exists()


# ── filter steps ──────────────────────────────────────────────────────────────

def test_filter_steps_first_is_all_cells(qc_adata):
    steps = _compute_filter_steps(qc_adata, THRESHOLDS)
    assert steps[0]["label"] == "All cells"
    assert steps[0]["n_remaining"] == qc_adata.n_obs
    assert steps[0]["n_removed"] == 0


def test_filter_steps_monotonically_decreasing(qc_adata):
    steps = _compute_filter_steps(qc_adata, THRESHOLDS)
    remaining = [s["n_remaining"] for s in steps]
    assert all(remaining[i] >= remaining[i + 1] for i in range(len(remaining) - 1))


def test_filter_steps_pct_removed_non_negative(qc_adata):
    steps = _compute_filter_steps(qc_adata, THRESHOLDS)
    assert all(s["pct_removed"] >= 0 for s in steps)


def test_filter_steps_empty_thresholds(qc_adata):
    """With empty thresholds, all cells should survive every step."""
    steps = _compute_filter_steps(qc_adata, {})
    assert all(s["n_remaining"] == qc_adata.n_obs for s in steps)


# ── metric stats ──────────────────────────────────────────────────────────────

def test_metric_stats_returns_four_metrics(qc_adata):
    stats = _compute_metric_stats(qc_adata)
    assert len(stats) == 4


def test_metric_stats_median_in_range(qc_adata):
    stats = _compute_metric_stats(qc_adata)
    for row in stats:
        assert row["min"] <= row["median"] <= row["max"]
        assert row["p5"] <= row["p95"]


# ── plot functions ────────────────────────────────────────────────────────────

def test_plot_violin_returns_figure(qc_adata):
    import matplotlib.pyplot as plt
    fig = plot_violin_with_thresholds(qc_adata, THRESHOLDS)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_scatter_returns_figure(qc_adata):
    import matplotlib.pyplot as plt
    fig = plot_scatter_with_thresholds(qc_adata, THRESHOLDS)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_histograms_returns_figure(qc_adata):
    import matplotlib.pyplot as plt
    fig = plot_qc_histograms(qc_adata, THRESHOLDS)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_knee_returns_figure(qc_adata):
    import matplotlib.pyplot as plt
    fig = plot_knee_curve(qc_adata)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_waterfall_returns_figure(qc_adata):
    import matplotlib.pyplot as plt
    steps = _compute_filter_steps(qc_adata, THRESHOLDS)
    fig = plot_filter_waterfall(steps)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_qc_saves_files(tmp_path, qc_adata):
    """Backward-compatible plot_qc() must still write PNGs."""
    from scanpy_workflow.qc.plots import plot_qc
    plot_qc(qc_adata, str(tmp_path))
    assert (tmp_path / "qc_violin.png").exists()
    assert (tmp_path / "qc_scatter_counts_vs_genes.png").exists()
