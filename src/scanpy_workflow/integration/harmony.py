import anndata as ad
import harmonypy as hm
from .evaluate import evaluate_batch_correction


def batch_correct_harmony(
    adata: ad.AnnData,
    batch_key: str = "sample",
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply Harmony batch correction on X_pca.

    Writes:
      - adata.obsm["X_pca_harmony"]
      - adata.uns["neighbors_use_rep"] = "X_pca_harmony"
    """
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca not found in adata.obsm. Run PCA first.")
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
    # harmonypy >= 0.2.0 returns Z_corr as (n_cells, n_pcs);
    # older versions returned (n_pcs, n_cells) requiring .T
    z = ho.Z_corr
    if z.shape[0] != adata.n_obs:
        z = z.T
    adata.obsm["X_pca_harmony"] = z
    adata.uns["neighbors_use_rep"] = "X_pca_harmony"

    if run_evaluation and output_dir:
        evaluate_batch_correction(adata, batch_key=batch_key, use_rep="X_pca_harmony",
                                   output_dir=output_dir)
    return adata
