import anndata as ad
import scanpy as sc


def scale(adata: ad.AnnData, max_value: float = 10.0) -> ad.AnnData:
    """
    Z-score scale gene expression. Sets adata.X to scaled data.
    Run after normalization and HVG selection, before PCA (non-scVI path).
    """
    sc.pp.scale(adata, max_value=max_value)
    return adata
