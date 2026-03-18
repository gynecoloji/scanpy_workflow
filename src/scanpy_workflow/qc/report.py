"""
Per-sample QC HTML report.

Usage
-----
    from scanpy_workflow.qc.report import generate_qc_report

    generate_qc_report(
        adata,
        sample_name  = "sample_A",
        thresholds   = {"min_genes": 200, "max_genes": 6000, ...},
        output_path  = "results/qc/sample_A/qc_report.html",
    )

The report is a single self-contained HTML file (plots embedded as base64 PNG).
"""
from __future__ import annotations

import base64
import io
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
from jinja2 import Environment, FileSystemLoader

from scanpy_workflow.qc.plots import (
    plot_violin_with_thresholds,
    plot_scatter_with_thresholds,
    plot_qc_histograms,
    plot_knee_curve,
    plot_filter_waterfall,
)


# ── public API ────────────────────────────────────────────────────────────────

def generate_qc_report(
    adata: ad.AnnData,
    sample_name: str,
    thresholds: dict,
    output_path: str,
) -> None:
    """Render a per-sample QC HTML report.

    Parameters
    ----------
    adata:
        AnnData **after** ``compute_qc_metrics()`` but **before** filtering,
        so the full cell distribution is visible.
    sample_name:
        Label used in the report header.
    thresholds:
        Dict of filter cutoffs, e.g.::

            {"min_genes": 200, "max_genes": 6000,
             "min_counts": 500, "max_counts": 30000,
             "max_pct_mito": 20.0, "max_pct_ribo": 50.0}

    output_path:
        Path to write the ``qc_report.html`` file.
    """
    required = ["n_genes_by_counts", "total_counts",
                "pct_counts_mito", "pct_counts_ribo"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Missing QC columns: {missing}. Run compute_qc_metrics() first."
        )

    filter_steps = _compute_filter_steps(adata, thresholds)
    metric_stats = _compute_metric_stats(adata)
    plots        = _render_plots(adata, thresholds, filter_steps)

    template_dir = Path(__file__).parent / "templates"
    env      = Environment(loader=FileSystemLoader(str(template_dir)),
                           autoescape=False)
    template = env.get_template("qc_report.html.j2")

    n_pass = filter_steps[-1]["n_remaining"] if filter_steps else adata.n_obs
    html   = template.render(
        sample_name  = sample_name,
        timestamp    = datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_cells_total= adata.n_obs,
        n_genes_total= adata.n_vars,
        n_cells_pass = n_pass,
        pct_pass     = round(100 * n_pass / max(adata.n_obs, 1), 1),
        thresholds   = thresholds,
        filter_steps = filter_steps,
        metric_stats = metric_stats,
        plots        = plots,
    )

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")


# ── helpers ───────────────────────────────────────────────────────────────────

def _compute_filter_steps(adata: ad.AnnData, thresholds: dict) -> list:
    """Return a list of dicts describing each filter step (sequential)."""
    obs = adata.obs

    steps_spec = [
        ("All cells",       None),
        ("min_genes",       lambda df: df["n_genes_by_counts"] >= thresholds.get("min_genes", 0)),
        ("max_genes",       lambda df: df["n_genes_by_counts"] <= thresholds.get("max_genes", 1e9)),
        ("min_counts",      lambda df: df["total_counts"]       >= thresholds.get("min_counts", 0)),
        ("max_counts",      lambda df: df["total_counts"]       <= thresholds.get("max_counts", 1e9)),
        ("max_pct_mito",    lambda df: df["pct_counts_mito"]    <= thresholds.get("max_pct_mito", 100)),
        ("max_pct_ribo",    lambda df: df["pct_counts_ribo"]    <= thresholds.get("max_pct_ribo", 100)),
    ]

    steps  = []
    mask   = pd.Series(True, index=obs.index)
    n_prev = adata.n_obs

    for label, condition in steps_spec:
        if condition is not None:
            if label not in thresholds and label not in (
                "min_genes", "max_genes", "min_counts",
                "max_counts", "max_pct_mito", "max_pct_ribo"
            ):
                continue
            mask = mask & condition(obs)

        n_now     = int(mask.sum())
        n_removed = n_prev - n_now
        pct_rem   = round(100 * n_removed / max(adata.n_obs, 1), 2)

        steps.append({
            "label":       label,
            "threshold":   thresholds.get(label, "—"),
            "n_remaining": n_now,
            "n_removed":   n_removed,
            "pct_removed": pct_rem,
        })
        n_prev = n_now

    return steps


def _compute_metric_stats(adata: ad.AnnData) -> list:
    """Per-metric summary statistics (median, mean, min, max, p5, p95)."""
    metrics = [
        ("n_genes_by_counts", "Genes per cell"),
        ("total_counts",       "UMI counts"),
        ("pct_counts_mito",    "% Mitochondrial"),
        ("pct_counts_ribo",    "% Ribosomal"),
    ]
    rows = []
    for col, label in metrics:
        if col not in adata.obs.columns:
            continue
        vals = adata.obs[col].dropna().values
        rows.append({
            "metric": label,
            "median": round(float(np.median(vals)), 2),
            "mean":   round(float(np.mean(vals)), 2),
            "min":    round(float(np.min(vals)), 2),
            "max":    round(float(np.max(vals)), 2),
            "p5":     round(float(np.percentile(vals, 5)), 2),
            "p95":    round(float(np.percentile(vals, 95)), 2),
        })
    return rows


def _fig_to_b64(fig: plt.Figure) -> str:
    """Encode a matplotlib Figure as a base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _render_plots(
    adata: ad.AnnData,
    thresholds: dict,
    filter_steps: list,
) -> dict:
    return {
        "violin":     _fig_to_b64(plot_violin_with_thresholds(adata, thresholds)),
        "scatter":    _fig_to_b64(plot_scatter_with_thresholds(adata, thresholds)),
        "histograms": _fig_to_b64(plot_qc_histograms(adata, thresholds)),
        "knee":       _fig_to_b64(plot_knee_curve(adata)),
        "waterfall":  _fig_to_b64(plot_filter_waterfall(filter_steps)),
    }
