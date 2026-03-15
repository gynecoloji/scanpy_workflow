import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional


def run_de(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 100,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Run differential expression analysis.

    Outputs a normalized DataFrame with columns:
      group, gene, score, logfoldchange, pval, pval_adj

    For logreg: pval and pval_adj are NaN (logreg returns coefficients only).
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs.")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    # Use only the groups that have results (logreg omits the reference group)
    # Use only the groups that have results (logreg omits the reference group)
    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    dfs = []
    for group in groups:
        result = sc.get.rank_genes_groups_df(adata, group=str(group))
        result = result.rename(columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "logfoldchange",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        })
        result["group"] = str(group)
        # Add missing columns (logreg does not return logfoldchange/pval/pval_adj)
        for col in ["logfoldchange", "pval", "pval_adj"]:
            if col not in result.columns:
                result[col] = np.nan
        dfs.append(result)

    df = pd.concat(dfs, ignore_index=True)
    df = df[["group", "gene", "score", "logfoldchange", "pval", "pval_adj"]]

    if method == "logreg":
        df["pval"] = np.nan
        df["pval_adj"] = np.nan

    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / "de_results.csv", index=False)

    return df
