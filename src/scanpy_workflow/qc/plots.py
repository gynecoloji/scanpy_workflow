from pathlib import Path
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt


def plot_qc(adata: ad.AnnData, output_dir: str) -> None:
    """Generate QC violin and scatter plots, save as PNG."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mito"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig(out / "qc_violin.png", bbox_inches="tight", dpi=100)
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=3,
        alpha=0.5,
    )
    ax.set_xlabel("Total counts")
    ax.set_ylabel("Genes by counts")
    fig.savefig(out / "qc_scatter_counts_vs_genes.png", bbox_inches="tight", dpi=100)
    plt.close(fig)
