from typing import List
import anndata as ad


def merge_samples(
    adatas: List[ad.AnnData],
    sample_names: List[str],
    batch_key: str = "sample",
) -> ad.AnnData:
    """
    Merge per-sample AnnData objects into a single AnnData.

    Creates adata.obs[batch_key] with the sample name for each cell,
    using anndata.concat with label=batch_key, keys=sample_names.
    """
    if len(adatas) != len(sample_names):
        raise ValueError(f"len(adatas)={len(adatas)} != len(sample_names)={len(sample_names)}")

    # Validate that all samples share the same var_names
    ref_var = adatas[0].var_names
    for i, a in enumerate(adatas[1:], start=1):
        if not ref_var.equals(a.var_names):
            raise ValueError(
                f"var_names of sample {i} do not match sample 0. "
                "All samples must have identical gene sets before merging."
            )

    merged = ad.concat(
        adatas,
        label=batch_key,
        keys=sample_names,
        merge="first",
        uns_merge="first",
    )
    return merged
