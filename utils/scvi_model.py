"""
Shared scVI VAE model training and inference utilities.
Used by both integration/scvi.py (batch correction) and imputation/scvi.py.
"""

import anndata as ad
from typing import Optional


def train_vae(
    adata: ad.AnnData,
    batch_key: str,
    layer: str = "counts",
    max_epochs: int = 400,
    early_stopping: bool = True,
    seed: int = 42,
):
    """
    Set up and train an scVI VAE model.

    Returns the trained SCVI model object.
    Input layer should be raw counts (layers['counts']).
    scVI's internal likelihood model requires raw counts.
    """
    import scvi
    scvi.settings.seed = seed

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
    )
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=max_epochs, early_stopping=early_stopping)
    return model
