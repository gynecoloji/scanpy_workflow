import anndata as ad
import scanpy as sc
import scipy.sparse
import numpy as np


def normalize(
    adata: ad.AnnData,
    method: str = "library_size",
    target_sum: float = 1e4,
) -> ad.AnnData:
    """
    Normalize raw counts and store three layers.

    Layers written:
      - layers["counts"]   : raw integer counts (copy of input X)
      - layers["norm"]     : library-size normalized counts (pre-log)
      - layers["log_norm"] : log1p of normalized counts

    adata.X is set to layers["log_norm"] on return.

    Args:
        adata: AnnData with raw counts in X.
        method: "library_size" (default) or "scran" (raises — use scran.R).
        target_sum: total counts per cell after normalization (library_size only).
    """
    if method == "scran":
        raise ValueError(
            "scran normalization must be run via scran.R. "
            "Use the R-backed normalize step in the Nextflow pipeline."
        )
    if method != "library_size":
        raise ValueError(f"Unknown normalization method: '{method}'. Use 'library_size' or 'scran'.")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Library-size normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.layers["norm"] = adata.X.copy()

    # Log-transform
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()
    # adata.X is already log_norm (set by log1p in place)

    return adata
