import anndata as ad
import scanpy as sc


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Compute QC metrics and add to adata.obs.

    Adds: n_genes_by_counts, total_counts, pct_counts_mito, pct_counts_ribo.
    Infers mt/ribo gene flags from var_names if not already in adata.var.
    """
    if "mt" not in adata.var.columns:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if "ribo" not in adata.var.columns:
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # Rename mito column for consistent naming (scanpy outputs pct_counts_mt, we use pct_counts_mito)
    adata.obs["pct_counts_mito"] = adata.obs["pct_counts_mt"]
    # pct_counts_ribo is already named correctly by scanpy (qc_vars=["ribo"] → pct_counts_ribo)
    return adata
