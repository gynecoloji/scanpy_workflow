import anndata as ad
from pathlib import Path
from scanpy_workflow.utils.scvi_model import train_vae


def batch_correct_scvi(
    adata: ad.AnnData,
    batch_key: str = "sample",
    max_epochs: int = 400,
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Apply scVI batch correction.

    Uses layers['counts'] (raw counts) as input.
    Writes:
      - adata.obsm["X_scVI"]         : low-dim latent embedding
      - adata.layers["scvi_denoised"] : denoised expression
      - adata.uns["neighbors_use_rep"] = "X_scVI"
    """
    model = train_vae(adata, batch_key=batch_key, layer="counts")

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_denoised"] = model.get_normalized_expression(library_size=1e4)
    adata.uns["neighbors_use_rep"] = "X_scVI"

    if run_evaluation and output_dir:
        from .evaluate import evaluate_batch_correction
        evaluate_batch_correction(adata, batch_key=batch_key,
                                   use_rep="X_scVI", output_dir=output_dir)
    return adata
