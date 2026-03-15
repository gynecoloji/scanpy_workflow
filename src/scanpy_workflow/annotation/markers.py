import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path
from typing import Optional


def rank_marker_genes(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 50,
    output_dir: Optional[str] = None,
) -> ad.AnnData:
    """
    Identify marker genes per cluster using rank_genes_groups.
    If output_dir given, writes marker_genes.csv.
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        groups = adata.obs[groupby].unique().tolist()
        dfs = [
            sc.get.rank_genes_groups_df(adata, group=g).assign(group=g)
            for g in groups
        ]
        pd.concat(dfs).to_csv(out / "marker_genes.csv", index=False)

    return adata
