import anndata as ad
from pathlib import Path
from scanpy_workflow.utils.scvi_model import train_vae


def impute_scvi(
    adata: ad.AnnData,
    batch_key: str = "sample",
    max_epochs: int = 400,
    model_dir: str = None,
    run_evaluation: bool = False,
    output_dir: str = None,
) -> ad.AnnData:
    """
    Impute using scVI (trains its own model on layers['counts'] (raw counts, as required by scVI)).

    Only valid when batch_correction.method != 'scvi' (enforced by config.py).
    Writes:
      - adata.layers["scvi_denoised"]
      - saves model to model_dir if provided
    """
    model = train_vae(adata, batch_key=batch_key, layer="counts", max_epochs=max_epochs)

    adata.layers["scvi_denoised"] = model.get_normalized_expression(library_size=1e4)

    if model_dir:
        Path(model_dir).mkdir(parents=True, exist_ok=True)
        model.save(model_dir, overwrite=True)

    if run_evaluation and output_dir:
        from .evaluate import evaluate_imputation
        import scvi as _scvi
        def _scvi_impute_fn(a: ad.AnnData) -> ad.AnnData:
            _scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key=batch_key)
            m = _scvi.model.SCVI(a)
            m.train(max_epochs=max_epochs, early_stopping=True)
            a.layers["scvi_denoised"] = m.get_normalized_expression(library_size=1e4)
            return a
        evaluate_imputation(adata, impute_fn=_scvi_impute_fn, output_dir=output_dir)
    return adata
