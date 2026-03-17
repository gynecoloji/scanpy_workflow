"""
Pathway scoring for single-cell data.

Three methods:
  aucell       — rank-based AUCell (via decoupler); robust to zero-inflation.
  ulm          — univariate linear model (via decoupler); requires weights.
  scanpy_score — scanpy.tl.score_genes; no extra dependency.

Results are stored in adata.obsm["pathway_scores"] (pd.DataFrame)
and written to:
  <output_dir>/pathway_scores.csv           (cells × pathways)
  <output_dir>/cluster_pathway_scores.csv   (clusters × pathways, mean)
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def run_pathway_scoring(
    adata: ad.AnnData,
    net: pd.DataFrame,
    method: str = "aucell",
    groupby: str = "leiden",
    min_n: int = 5,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Score cells against gene sets and return a cells × pathways DataFrame.

    Parameters
    ----------
    adata:
        AnnData object; adata.X should be log-normalised counts.
    net:
        Tidy gene-set DataFrame with at minimum columns ``source`` and ``target``.
        A ``weight`` column is used by the ulm method; ignored by aucell/scanpy_score.
    method:
        "aucell" | "ulm" | "scanpy_score"
    groupby:
        obs column used to compute per-cluster mean summary table.
    min_n:
        Minimum number of genes from a gene set that must be present in adata.var_names.
        Gene sets with fewer overlapping genes are dropped.
    output_dir:
        Directory to write CSV outputs; created if absent.
    """
    if "source" not in net.columns or "target" not in net.columns:
        raise ValueError("net must have 'source' and 'target' columns.")

    if method in {"aucell", "ulm"}:
        scores = _score_decoupler(adata, net, method, min_n)
    elif method == "scanpy_score":
        scores = _score_scanpy(adata, net, min_n)
    else:
        raise ValueError(
            f"Unknown method: '{method}'. Valid: aucell, ulm, scanpy_score"
        )

    if scores.empty:
        raise RuntimeError(
            "No gene sets survived the min_n filter — "
            "check that gene symbols match adata.var_names."
        )

    adata.obsm["pathway_scores"] = scores

    if output_dir:
        _write_outputs(scores, adata, groupby, output_dir)

    return scores


# ── scoring back-ends ──────────────────────────────────────────────────────────

def _score_decoupler(
    adata: ad.AnnData,
    net: pd.DataFrame,
    method: str,
    min_n: int,
) -> pd.DataFrame:
    try:
        import decoupler as dc
    except ImportError:
        raise ImportError(
            f"decoupler-py is required for method='{method}'. "
            "Install with: pip install decoupler-py"
        )

    net = net.copy()
    # decoupler expects 'source' and 'target' columns; ensure they exist
    weight_col = "weight" if "weight" in net.columns else None

    if method == "aucell":
        dc.run_aucell(
            mat=adata,
            net=net,
            source="source",
            target="target",
            min_n=min_n,
            use_raw=False,
        )
        scores: pd.DataFrame = adata.obsm["aucell_estimate"]

    elif method == "ulm":
        dc.run_ulm(
            mat=adata,
            net=net,
            source="source",
            target="target",
            weight=weight_col,
            min_n=min_n,
            use_raw=False,
        )
        scores = adata.obsm["ulm_estimate"]

    # Clean intermediate obsm keys so only pathway_scores remains
    for key in ("aucell_estimate", "aucell_pvals", "ulm_estimate", "ulm_pvals"):
        adata.obsm.pop(key, None)

    return scores


def _score_scanpy(
    adata: ad.AnnData,
    net: pd.DataFrame,
    min_n: int,
) -> pd.DataFrame:
    """Use scanpy.tl.score_genes (control-gene subtracted mean)."""
    gene_sets: dict[str, list[str]] = (
        net.groupby("source")["target"]
        .apply(list)
        .to_dict()
    )

    score_cols: list[str] = []
    for gs_name, genes in gene_sets.items():
        available = [g for g in genes if g in adata.var_names]
        if len(available) < min_n:
            continue
        # Sanitize name for use as obs column key (max 64 chars)
        key = re.sub(r"[^A-Za-z0-9_]", "_", gs_name)[:64]
        # Guarantee uniqueness in the unlikely event of collisions
        if key in score_cols:
            key = f"{key}_{len(score_cols)}"
        sc.tl.score_genes(adata, gene_list=available, score_name=key)
        score_cols.append(key)

    if not score_cols:
        return pd.DataFrame(index=adata.obs_names)

    scores = adata.obs[score_cols].copy()
    # Remove from obs — scores live in obsm["pathway_scores"]
    adata.obs.drop(columns=score_cols, inplace=True)
    scores.columns = [c.replace("_", " ").title() if c not in score_cols else c
                      for c in scores.columns]
    # Restore original column names (avoid double rename)
    scores.columns = score_cols
    return scores


# ── output helpers ─────────────────────────────────────────────────────────────

def _write_outputs(
    scores: pd.DataFrame,
    adata: ad.AnnData,
    groupby: str,
    output_dir: str,
) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    scores.to_csv(out / "pathway_scores.csv")

    # Per-cluster mean summary
    if groupby in adata.obs.columns:
        cluster_labels = adata.obs[groupby].values
        summary = (
            scores.assign(_cluster=cluster_labels)
            .groupby("_cluster")
            .mean()
        )
        summary.index.name = groupby
        summary.to_csv(out / "cluster_pathway_scores.csv")
