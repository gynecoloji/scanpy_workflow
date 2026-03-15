import anndata as ad
import scanpy as sc


def filter_cells_genes(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    min_counts: int = 500,
    max_counts: int = 30000,
    max_pct_mito: float = 20.0,
    max_pct_ribo: float = 50.0,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Filter cells and genes based on QC cutoffs.
    Requires compute_qc_metrics() to have been run first.
    """
    required = ["n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ribo"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing QC columns: {missing}. Run compute_qc_metrics() first.")

    cell_mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["pct_counts_mito"] <= max_pct_mito)
        & (adata.obs["pct_counts_ribo"] <= max_pct_ribo)
    )
    adata = adata[cell_mask].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata
