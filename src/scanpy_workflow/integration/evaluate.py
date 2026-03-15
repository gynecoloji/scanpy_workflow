"""Batch correction evaluation utilities."""
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


def evaluate_batch_correction(
    adata: ad.AnnData,
    batch_key: str,
    use_rep: str,
    output_dir: str,
    label_key: str = None,
) -> dict:
    """
    Evaluate batch correction using scib-metrics.

    Computes: kBET approximation, iLISI, cLISI, ASW batch.
    Writes UMAP post-correction plot and metrics JSON to output_dir.
    Returns dict of metric scores.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    scores = {}

    try:
        import scib_metrics
        from scib_metrics.benchmark import Benchmarker
        bm_kwargs = dict(
            adata=adata,
            batch_key=batch_key,
            embedding_obsm_keys=[use_rep],
        )
        if label_key is not None:
            bm_kwargs["label_key"] = label_key
        try:
            bm = Benchmarker(**bm_kwargs)
            bm.benchmark()
            results = bm.get_results(min_max_scale=False)
            scores = results.to_dict()
        except Exception as e:
            scores["error"] = str(e)
    except ImportError as e:
        scores["error"] = str(e)

    # UMAP comparison plot
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=batch_key, show=False)
    plt.savefig(out / "umap_post_correction.png", bbox_inches="tight", dpi=100)
    plt.close()

    import json
    with open(out / "batch_correction_metrics.json", "w") as f:
        json.dump(scores, f, indent=2, default=str)

    return scores
