"""Batch correction evaluation utilities."""
import anndata as ad
from pathlib import Path


def evaluate_batch_correction(
    adata: ad.AnnData,
    batch_key: str = "sample",
    use_rep: str = "X_pca",
    output_dir: str = None,
) -> None:
    """
    Evaluate batch correction performance.

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object with batch correction applied.
    batch_key : str
        Column name in adata.obs containing batch labels. Default: "sample".
    use_rep : str
        Representation to use for evaluation. Default: "X_pca".
    output_dir : str, optional
        Output directory for evaluation results.
    """
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    # TODO: Implement batch correction evaluation metrics
