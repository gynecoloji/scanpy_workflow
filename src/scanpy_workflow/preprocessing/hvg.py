import anndata as ad
import scanpy as sc


def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 3000,
    flavor: str = "seurat_v3",
) -> ad.AnnData:
    """
    Select highly variable genes and add 'highly_variable' flag to adata.var.

    For seurat_v3: requires layers["counts"] (raw counts).
    For seurat / cell_ranger: uses adata.X (log-normalized).
    """
    if flavor == "seurat_v3":
        if "counts" not in adata.layers:
            raise ValueError(
                "layers['counts'] required for seurat_v3 flavor. "
                "Run normalize() first to store raw counts."
            )
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            layer="counts",
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
        )
    return adata
