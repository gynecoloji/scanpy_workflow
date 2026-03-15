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

    merged = ad.concat(
        adatas,
        label=batch_key,
        keys=sample_names,
        merge="first",
        uns_merge="first",
    )
    return merged
