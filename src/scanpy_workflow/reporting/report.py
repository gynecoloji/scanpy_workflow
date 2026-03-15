import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import json
from datetime import datetime
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from typing import Optional


def generate_report(
    adata: ad.AnnData,
    output_dir: str,
    params: dict,
    cluster_key: str = "louvain",
    cluster_metrics: Optional[dict] = None,
) -> None:
    """Generate an HTML summary report."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # UMAP plot
    umap_plot = None
    if "X_umap" in adata.obsm and cluster_key in adata.obs.columns:
        sc.pl.umap(adata, color=cluster_key, show=False)
        umap_file = out / "umap_clusters.png"
        plt.savefig(umap_file, bbox_inches="tight", dpi=100)
        plt.close()
        umap_plot = "umap_clusters.png"

    samples = list(adata.obs.get("sample", adata.obs.get("batch", pd.Series(["unknown"]))).unique())
    n_clusters = adata.obs[cluster_key].nunique() if cluster_key in adata.obs.columns else "N/A"

    template_dir = Path(__file__).parent / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)))
    template = env.get_template("report.html.j2")

    html = template.render(
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        samples=params.get("samples", samples),
        n_clusters=n_clusters,
        umap_plot=umap_plot,
        cluster_metrics=cluster_metrics or {},
        params_json=json.dumps(params, indent=2, default=str),
    )
    (out / "report.html").write_text(html)
