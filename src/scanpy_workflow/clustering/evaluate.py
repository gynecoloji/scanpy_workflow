import anndata as ad
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
from pathlib import Path

SUBSAMPLE_N = 20_000
RANDOM_STATE = 42


def evaluate_clustering(
    adata: ad.AnnData,
    cluster_key: str,
    output_dir: str,
) -> dict:
    """
    Compute clustering evaluation metrics.

    Silhouette and Davies-Bouldin: subsampled to <= SUBSAMPLE_N cells.
    Calinski-Harabasz: full dataset.
    """
    from sklearn.metrics import (
        silhouette_score,
        davies_bouldin_score,
        calinski_harabasz_score,
    )
    import scanpy as sc

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    use_rep = adata.uns.get("neighbors_use_rep", "X_pca")
    embedding = adata.obsm.get(use_rep)
    if embedding is None:
        embedding = adata.obsm.get("X_pca")
    if embedding is None:
        raise ValueError("No valid embedding found for clustering evaluation.")

    labels = adata.obs[cluster_key].astype(str).values

    # Subsample for O(n²) metrics
    n = adata.n_obs
    if n > SUBSAMPLE_N:
        rng = np.random.default_rng(RANDOM_STATE)
        idx = rng.choice(n, size=SUBSAMPLE_N, replace=False)
        emb_sub = embedding[idx]
        lab_sub = labels[idx]
    else:
        emb_sub = embedding
        lab_sub = labels

    scores = {}
    if len(np.unique(lab_sub)) > 1:
        scores["silhouette"] = float(silhouette_score(emb_sub, lab_sub))
        scores["davies_bouldin"] = float(davies_bouldin_score(emb_sub, lab_sub))
        scores["calinski_harabasz"] = float(calinski_harabasz_score(embedding, labels))
    else:
        scores = {"silhouette": None, "davies_bouldin": None, "calinski_harabasz": None,
                  "note": "Only one cluster found; metrics not computed."}

    with open(out / "clustering_metrics.json", "w") as f:
        json.dump(scores, f, indent=2)

    # Marker dot plot
    try:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=5)
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        plt.savefig(out / "marker_dotplot.png", bbox_inches="tight", dpi=100)
        plt.close()
    except Exception:
        pass  # Skip if not enough cells per cluster

    return scores
