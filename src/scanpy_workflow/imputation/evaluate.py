import anndata as ad
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse
from pathlib import Path


def evaluate_imputation(
    adata: ad.AnnData,
    impute_fn: callable,
    output_dir: str,
    mask_fraction: float = 0.10,
    random_state: int = 42,
) -> dict:
    """
    In-silico dropout evaluation for imputation.

    Makes a deep copy of adata, masks mask_fraction of non-zero values in
    layers["norm"] of the copy, calls impute_fn(adata_copy) to get imputed
    results on the masked data, then compares imputed values vs original true
    values at masked positions.

    The production adata is never modified.
    Returns dict with pearson_r, spearman_r.
    """
    import copy
    from scipy.stats import pearsonr, spearmanr
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Make a deep copy to avoid modifying production adata
    adata_copy = copy.deepcopy(adata)

    # Get the norm layer from the copy for masking
    layer = adata_copy.layers["norm"]
    if scipy.sparse.issparse(layer):
        layer_dense = layer.toarray().copy()
    else:
        layer_dense = np.asarray(layer).copy()

    # Record original (true) values before masking
    rng = np.random.default_rng(random_state)
    nonzero_i, nonzero_j = np.nonzero(layer_dense)
    n_mask = int(len(nonzero_i) * mask_fraction)
    idx = rng.choice(len(nonzero_i), size=n_mask, replace=False)
    mask_i, mask_j = nonzero_i[idx], nonzero_j[idx]

    true_values = layer_dense[mask_i, mask_j].copy()

    # Apply mask to the copy's norm layer
    masked_layer = layer_dense.copy()
    masked_layer[mask_i, mask_j] = 0.0
    adata_copy.layers["norm"] = masked_layer

    # Run imputation on the masked copy
    imputed_adata = impute_fn(adata_copy)

    # Get imputed values at masked positions
    after_layer = imputed_adata.layers.get("imputed", imputed_adata.X)
    if scipy.sparse.issparse(after_layer):
        after_dense = after_layer.toarray()
    else:
        after_dense = np.asarray(after_layer)

    imputed_values = after_dense[mask_i, mask_j]

    pearson_r, _ = pearsonr(true_values, imputed_values)
    spearman_r, _ = spearmanr(true_values, imputed_values)
    scores = {"pearson_r": float(pearson_r), "spearman_r": float(spearman_r)}

    import json
    with open(out / "imputation_metrics.json", "w") as f:
        json.dump(scores, f, indent=2)

    # Distribution plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(layer_dense.flatten(), bins=50, alpha=0.7, label="before")
    axes[0].set_title("Before imputation")
    axes[1].hist(after_dense.flatten(), bins=50, alpha=0.7, label="after", color="orange")
    axes[1].set_title("After imputation")
    fig.savefig(out / "imputation_distribution.png", bbox_inches="tight", dpi=100)
    plt.close(fig)

    return scores
