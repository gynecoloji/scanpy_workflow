import anndata as ad
import bbknn


def batch_correct_bbknn(
    adata: ad.AnnData,
    batch_key: str = "sample",
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply BBKNN batch correction.

    Input: adata with X_pca in obsm.
    Output:
      - adata.obsp["connectivities"] and adata.obsp["distances"] written directly.
      - adata.uns["neighbors_use_rep"] = "skip" (pp.neighbors must be skipped).

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object with X_pca in obsm and batch_key in obs.
    batch_key : str
        Column name in adata.obs containing batch labels. Default: "sample".
    run_evaluation : bool
        Whether to run batch correction evaluation. Default: False.
    output_dir : str, optional
        Output directory for evaluation results if run_evaluation is True.

    Returns
    -------
    ad.AnnData
        Modified adata with BBKNN batch correction applied.

    Raises
    ------
    ValueError
        If X_pca is not found in adata.obsm.
        If batch_key is not found in adata.obs.
    """
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca not found in adata.obsm. Run PCA first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    bbknn.bbknn(adata, batch_key=batch_key, use_rep="X_pca")
    adata.uns["neighbors_use_rep"] = "skip"

    if run_evaluation and output_dir:
        from .evaluate import evaluate_batch_correction

        evaluate_batch_correction(
            adata, batch_key=batch_key, use_rep="X_pca", output_dir=output_dir
        )
    return adata
