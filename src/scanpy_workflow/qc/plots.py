"""
QC plot utilities.

Basic interface (used by the QC Nextflow step):
    plot_qc(adata, output_dir)          – violin + scatter, saves PNGs

Report interface (used by the QC report generator, return Figure objects):
    plot_violin_with_thresholds(...)
    plot_scatter_with_thresholds(...)
    plot_qc_histograms(...)
    plot_knee_curve(...)
    plot_filter_waterfall(...)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad


# ── basic interface (backward-compatible) ─────────────────────────────────────

def plot_qc(adata: ad.AnnData, output_dir: str) -> None:
    """Generate QC violin and scatter plots, save as PNG."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    fig = plot_violin_with_thresholds(adata, thresholds={})
    fig.savefig(out / "qc_violin.png", bbox_inches="tight", dpi=100)
    plt.close(fig)

    fig = plot_scatter_with_thresholds(adata, thresholds={})
    fig.savefig(out / "qc_scatter_counts_vs_genes.png", bbox_inches="tight", dpi=100)
    plt.close(fig)


# ── report-quality plots (return Figure) ──────────────────────────────────────

def plot_violin_with_thresholds(
    adata: ad.AnnData,
    thresholds: dict,
) -> plt.Figure:
    """4-panel violin plots with dashed threshold lines."""
    panels = [
        ("n_genes_by_counts", "Genes per cell",
         [("min_genes", "tab:red"), ("max_genes", "tab:orange")]),
        ("total_counts", "UMI counts",
         [("min_counts", "tab:red"), ("max_counts", "tab:orange")]),
        ("pct_counts_mito", "% Mitochondrial",
         [("max_pct_mito", "tab:orange")]),
        ("pct_counts_ribo", "% Ribosomal",
         [("max_pct_ribo", "tab:orange")]),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for ax, (col, label, thresh_keys) in zip(axes, panels):
        if col not in adata.obs.columns:
            ax.set_visible(False)
            continue
        vals = adata.obs[col].dropna().values
        parts = ax.violinplot(vals, positions=[0], widths=0.7,
                               showmedians=True, showextrema=True)
        for pc in parts.get("bodies", []):
            pc.set_facecolor("steelblue")
            pc.set_alpha(0.6)
        ax.set_xticks([])
        ax.set_ylabel(label)
        ax.set_title(label, fontsize=9)

        for key, color in thresh_keys:
            val = thresholds.get(key)
            if val is not None:
                ax.axhline(val, color=color, linestyle="--",
                           linewidth=1.5, label=f"{key}={val}")
        if any(thresholds.get(k) is not None for k, _ in thresh_keys):
            ax.legend(fontsize=7)

    fig.suptitle("QC Metric Distributions", fontsize=12, y=1.01)
    fig.tight_layout()
    return fig


def plot_scatter_with_thresholds(
    adata: ad.AnnData,
    thresholds: dict,
) -> plt.Figure:
    """Scatter: total_counts vs n_genes, coloured by % mito. Threshold lines overlaid."""
    fig, ax = plt.subplots(figsize=(6, 5))

    x = adata.obs.get("total_counts",      pd.Series(np.zeros(adata.n_obs)))
    y = adata.obs.get("n_genes_by_counts", pd.Series(np.zeros(adata.n_obs)))
    c = adata.obs.get("pct_counts_mito",   pd.Series(np.zeros(adata.n_obs)))

    sc_plot = ax.scatter(x, y, c=c, cmap="Reds", s=3, alpha=0.5,
                         vmin=0, vmax=max(float(c.max()), 1))
    plt.colorbar(sc_plot, ax=ax, label="% Mito", shrink=0.8)

    _lines = {
        "min_counts": ("v", "tab:red"),
        "max_counts": ("v", "tab:orange"),
        "min_genes":  ("h", "tab:red"),
        "max_genes":  ("h", "tab:orange"),
    }
    for key, (orient, color) in _lines.items():
        val = thresholds.get(key)
        if val is not None:
            fn = ax.axvline if orient == "v" else ax.axhline
            fn(val, color=color, linestyle="--", linewidth=1.2, label=f"{key}={val}")

    if any(thresholds.get(k) is not None for k in _lines):
        ax.legend(fontsize=8)

    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Genes detected")
    ax.set_title("Counts vs Genes (colour = % mito)")
    fig.tight_layout()
    return fig


def plot_qc_histograms(
    adata: ad.AnnData,
    thresholds: dict,
) -> plt.Figure:
    """4-panel histograms with vertical threshold lines."""
    panels = [
        ("n_genes_by_counts", "Genes per cell",
         [("min_genes", "tab:red"), ("max_genes", "tab:orange")]),
        ("total_counts", "UMI counts",
         [("min_counts", "tab:red"), ("max_counts", "tab:orange")]),
        ("pct_counts_mito", "% Mitochondrial",
         [("max_pct_mito", "tab:orange")]),
        ("pct_counts_ribo", "% Ribosomal",
         [("max_pct_ribo", "tab:orange")]),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    for ax, (col, label, thresh_keys) in zip(axes, panels):
        if col not in adata.obs.columns:
            ax.set_visible(False)
            continue
        vals = adata.obs[col].dropna().values
        ax.hist(vals, bins=60, color="steelblue", alpha=0.75, edgecolor="none")
        ax.set_xlabel(label)
        ax.set_ylabel("Cells")
        ax.set_title(label, fontsize=9)

        for key, color in thresh_keys:
            val = thresholds.get(key)
            if val is not None:
                ax.axvline(val, color=color, linestyle="--",
                           linewidth=1.5, label=f"{key}={val}")
        if any(thresholds.get(k) is not None for k, _ in thresh_keys):
            ax.legend(fontsize=7)

    fig.suptitle("QC Metric Histograms with Filter Thresholds", fontsize=12, y=1.01)
    fig.tight_layout()
    return fig


def plot_knee_curve(adata: ad.AnnData) -> plt.Figure:
    """Cell rank vs total UMI counts (log–log). Helps identify empty droplets."""
    sorted_counts = np.sort(adata.obs["total_counts"].values)[::-1]
    rank = np.arange(1, len(sorted_counts) + 1)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(rank, sorted_counts, lw=1.5, color="steelblue")
    ax.set_xlabel("Cell rank (by total counts)")
    ax.set_ylabel("Total UMI counts")
    ax.set_title("Knee plot")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    return fig


def plot_filter_waterfall(filter_steps: list) -> plt.Figure:
    """Bar chart showing cells remaining after each sequential filter."""
    labels    = [s["label"] for s in filter_steps]
    remaining = [s["n_remaining"] for s in filter_steps]
    n_start   = remaining[0] if remaining else 1

    colors = ["#4C72B0"] + [
        "#55A868" if r == remaining[i] else "#C44E52"
        for i, r in enumerate(remaining[1:], 1)
    ]

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.4), 4))
    bars = ax.bar(labels, remaining, color=colors, edgecolor="white", width=0.6)

    for bar, step in zip(bars, filter_steps):
        n   = step["n_remaining"]
        pct = step.get("pct_removed", 0.0)
        txt = f"{n:,}" if pct == 0 else f"{n:,}\n(−{pct:.1f}%)"
        ax.text(bar.get_x() + bar.get_width() / 2,
                n + n_start * 0.01, txt,
                ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Cells remaining")
    ax.set_title("Cells remaining after each filter")
    ax.set_ylim(0, n_start * 1.18)
    plt.xticks(rotation=25, ha="right")
    fig.tight_layout()
    return fig
