import anndata as ad
import scrublet as scr
import numpy as np


def detect_doublets(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.05,
    filter_doublets: bool = False,
    random_state: int = 42,
) -> ad.AnnData:
    """
    Detect doublets with Scrublet.

    Adds adata.obs['doublet_score'] and adata.obs['predicted_doublet'].
    If filter_doublets=True, removes predicted doublets.
    Uses raw counts (layers['counts'] if available, else adata.X).
    """
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X

    scrub = scr.Scrublet(
        counts_matrix=counts,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state,
    )
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    if filter_doublets:
        adata = adata[~adata.obs["predicted_doublet"]].copy()

    return adata
