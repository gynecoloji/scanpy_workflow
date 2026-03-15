import anndata as ad
import scanpy as sc
import numpy as np
import scipy.sparse


def run_pca(
    adata: ad.AnnData,
    n_comps: int = 50,
    scvi_path: bool = False,
) -> ad.AnnData:
    """
    Run PCA and add adata.obsm['X_pca'].

    scvi_path=True: asserts adata.X == layers['log_norm'] (scaling was skipped).
    scvi_path=False: assumes adata.X is scaled data.
    """
    if scvi_path:
        if "log_norm" not in adata.layers:
            raise ValueError("layers['log_norm'] required for scVI PCA path.")
        x = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.asarray(adata.X)
        ln = adata.layers["log_norm"].toarray() if scipy.sparse.issparse(adata.layers["log_norm"]) else np.asarray(adata.layers["log_norm"])
        assert np.allclose(x, ln, atol=1e-5), (
            "adata.X must equal layers['log_norm'] on the scVI path. "
            "Ensure scale() was NOT called before PCA when using scVI batch correction."
        )

    sc.pp.pca(adata, n_comps=min(n_comps, adata.n_obs - 1, adata.n_vars - 1),
              use_highly_variable=True if "highly_variable" in adata.var.columns else False)
    return adata
