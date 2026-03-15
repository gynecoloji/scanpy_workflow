import anndata as ad
import scanpy as sc


def cluster(
    adata: ad.AnnData,
    algorithm: str = "leiden",
    resolution: float = 0.5,
    n_neighbors: int = 15,
) -> ad.AnnData:
    """
    Build neighbor graph (unless BBKNN already did it) and run Leiden/Louvain.

    Reads adata.uns['neighbors_use_rep']:
      - "skip": BBKNN wrote the graph; skip pp.neighbors.
      - anything else: call pp.neighbors(use_rep=...).
    """
    use_rep = adata.uns.get("neighbors_use_rep", "X_pca")

    if use_rep != "skip":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    if algorithm == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
    elif algorithm == "louvain":
        sc.tl.louvain(adata, resolution=resolution)
    else:
        raise ValueError(
            f"Unknown clustering algorithm: '{algorithm}'. Use 'leiden' or 'louvain'."
        )
    return adata
