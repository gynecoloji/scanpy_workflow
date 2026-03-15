import anndata as ad
import celltypist
from typing import Optional
from pathlib import Path


def annotate_celltypist(
    adata: ad.AnnData,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    output_dir: Optional[str] = None,
) -> ad.AnnData:
    """
    Automated cell type annotation with CellTypist.

    Adds adata.obs['celltypist_cell_type'] and adata.obs['celltypist_conf_score'].
    CellTypist expects log-normalized data in adata.X.
    """
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
    )
    adata = predictions.to_adata()
    adata.obs["celltypist_cell_type"] = adata.obs.get(
        "majority_voting", adata.obs.get("predicted_labels", "unknown")
    )
    adata.obs["celltypist_conf_score"] = adata.obs.get("conf_score", 0.0)

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        adata.obs[["celltypist_cell_type", "celltypist_conf_score"]].to_csv(
            out / "celltypist_predictions.csv"
        )
    return adata
