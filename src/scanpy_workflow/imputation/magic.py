import anndata as ad
import numpy as np
import scipy.sparse


def impute_magic(
    adata: ad.AnnData,
    solver: str = "exact",
    t: int = 3,
    run_evaluation: bool = False,
    output_dir: str = None,
    random_state: int = 42,
) -> ad.AnnData:
    """
    Apply MAGIC imputation.

    Operates on layers['norm'] (library-size normalized, pre-log) per
    MAGIC's recommended usage. Sets adata.X to imputed values.
    """
    import magic

    if "norm" not in adata.layers:
        raise ValueError("layers['norm'] required for MAGIC. Run normalize() first.")

    norm = adata.layers["norm"]
    if scipy.sparse.issparse(norm):
        norm = norm.toarray()

    magic_op = magic.MAGIC(solver=solver, t=t, random_state=random_state, verbose=0)
    imputed = magic_op.fit_transform(norm)

    adata.layers["imputed"] = imputed
    adata.X = scipy.sparse.csr_matrix(imputed) if scipy.sparse.issparse(adata.layers["norm"]) else imputed

    if run_evaluation and output_dir:
        from .evaluate import evaluate_imputation
        evaluate_imputation(adata,
                             impute_fn=lambda a: impute_magic(a, run_evaluation=False),
                             output_dir=output_dir)
    return adata
