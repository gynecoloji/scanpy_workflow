"""
PCA evaluation utilities.

Outputs written to output_dir:
  pca_scree.png          — per-PC explained variance (bar) + cumulative (line)
  pca_loadings.png       — top gene loadings for PC1–PC4
  pca_metrics.json       — n_comps_used, variance_explained (list),
                           cumvar_80pct, cumvar_90pct, elbow_pc
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad


def evaluate_pca(
    adata: ad.AnnData,
    output_dir: str,
) -> dict:
    """Compute PCA quality metrics and write plots + JSON to output_dir.

    Requires ``sc.pp.pca`` to have been run (``adata.uns['pca']`` present).

    Returns a dict with keys:
        n_comps_used, variance_explained, cumvar_80pct, cumvar_90pct, elbow_pc
    """
    if "pca" not in adata.uns or "variance_ratio" not in adata.uns["pca"]:
        raise ValueError(
            "PCA results not found in adata.uns['pca']. Run sc.pp.pca() first."
        )

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    var_ratio = np.array(adata.uns["pca"]["variance_ratio"])  # fraction per PC
    cumvar    = np.cumsum(var_ratio)
    n_comps   = len(var_ratio)

    metrics = {
        "n_comps_used":     n_comps,
        "variance_explained": [round(float(v), 6) for v in var_ratio],
        "cumvar_80pct":     int(_comps_for_threshold(cumvar, 0.80)),
        "cumvar_90pct":     int(_comps_for_threshold(cumvar, 0.90)),
        "elbow_pc":         int(_detect_elbow(var_ratio)),
    }

    with open(out / "pca_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    _plot_scree(var_ratio, cumvar, metrics, out)

    if "PCs" in adata.varm:
        _plot_loadings(adata, out)

    return metrics


# ── metrics helpers ────────────────────────────────────────────────────────────

def _comps_for_threshold(cumvar: np.ndarray, threshold: float) -> int:
    """Return 1-indexed number of PCs needed to reach `threshold` cumulative variance."""
    idx = np.searchsorted(cumvar, threshold)
    return min(int(idx) + 1, len(cumvar))


def _detect_elbow(var_ratio: np.ndarray) -> int:
    """Return 1-indexed elbow PC via maximum perpendicular distance (geometric method)."""
    n = len(var_ratio)
    if n < 3:
        return n
    x = np.linspace(0.0, 1.0, n)
    y_min, y_max = var_ratio.min(), var_ratio.max()
    y = (var_ratio - y_min) / (y_max - y_min + 1e-12)
    # Perpendicular distance from each point to the line connecting first and last
    dx = x[-1] - x[0]
    dy = y[-1] - y[0]
    length = np.hypot(dx, dy) + 1e-12
    dist = np.abs(dy * x - dx * y + x[-1] * y[0] - y[-1] * x[0]) / length
    return int(np.argmax(dist)) + 1  # 1-indexed


# ── plot helpers ───────────────────────────────────────────────────────────────

def _plot_scree(
    var_ratio: np.ndarray,
    cumvar: np.ndarray,
    metrics: dict,
    out: Path,
) -> None:
    n = len(var_ratio)
    pcs = np.arange(1, n + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))

    # ── left: per-PC bar chart ──
    ax1.bar(pcs, var_ratio * 100, color="steelblue", alpha=0.8, width=0.8)
    elbow = metrics["elbow_pc"]
    ax1.axvline(elbow, color="tab:red", linestyle="--", linewidth=1.5,
                label=f"Elbow PC{elbow}")
    ax1.set_xlabel("Principal component")
    ax1.set_ylabel("Variance explained (%)")
    ax1.set_title("Scree plot")
    ax1.legend(fontsize=9)
    if n > 30:
        ax1.set_xlim(0.5, min(n + 0.5, 60))  # show at most first 60 PCs

    # ── right: cumulative variance ──
    ax2.plot(pcs, cumvar * 100, color="steelblue", lw=2)
    ax2.fill_between(pcs, cumvar * 100, alpha=0.15, color="steelblue")
    for threshold, label, color in [
        (0.80, "80%", "tab:orange"),
        (0.90, "90%", "tab:red"),
    ]:
        n_needed = _comps_for_threshold(cumvar, threshold)
        ax2.axhline(threshold * 100, color=color, linestyle="--",
                    linewidth=1.2, label=f"{label} @ PC{n_needed}")
        ax2.axvline(n_needed, color=color, linestyle=":", linewidth=1.0)

    ax2.set_xlabel("Number of principal components")
    ax2.set_ylabel("Cumulative variance explained (%)")
    ax2.set_title("Cumulative explained variance")
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 101)

    fig.tight_layout()
    fig.savefig(out / "pca_scree.png", bbox_inches="tight", dpi=100)
    plt.close(fig)


def _plot_loadings(adata: ad.AnnData, out: Path, n_pcs: int = 4, top_n: int = 15) -> None:
    """Bar charts of top gene loadings for the first n_pcs PCs."""
    loadings = adata.varm["PCs"]          # genes × n_pcs
    gene_names = np.array(adata.var_names)
    n_pcs = min(n_pcs, loadings.shape[1])

    fig, axes = plt.subplots(1, n_pcs, figsize=(n_pcs * 4, 4), sharey=False)
    if n_pcs == 1:
        axes = [axes]

    for pc_idx, ax in enumerate(axes):
        loading = loadings[:, pc_idx]
        top_idx = np.argsort(np.abs(loading))[-top_n:][::-1]
        vals  = loading[top_idx]
        genes = gene_names[top_idx]
        colors = ["tab:blue" if v >= 0 else "tab:red" for v in vals]
        ax.barh(range(top_n), vals[::-1], color=colors[::-1], alpha=0.8)
        ax.set_yticks(range(top_n))
        ax.set_yticklabels(genes[::-1], fontsize=7)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(f"PC{pc_idx + 1}", fontsize=10)
        ax.set_xlabel("Loading")

    fig.suptitle(f"Top {top_n} gene loadings", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(out / "pca_loadings.png", bbox_inches="tight", dpi=100)
    plt.close(fig)
