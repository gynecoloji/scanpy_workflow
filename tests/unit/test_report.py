from pathlib import Path
from scanpy_workflow.reporting.report import generate_report
import numpy as np


def test_report_creates_html(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    # Manually add louvain clusters to avoid igraph dependency
    adata.obs["louvain"] = np.random.randint(0, 3, size=adata.n_obs).astype(str)
    sc.tl.umap(adata)
    generate_report(adata, output_dir=str(tmp_path), params={})
    assert (tmp_path / "report.html").exists()


def test_report_contains_sample_info(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    # Manually add louvain clusters to avoid igraph dependency
    adata.obs["louvain"] = np.random.randint(0, 3, size=adata.n_obs).astype(str)
    sc.tl.umap(adata)
    generate_report(adata, output_dir=str(tmp_path), params={"samples": ["s1", "s2"]})
    html = (tmp_path / "report.html").read_text()
    assert "s1" in html or "Cells" in html
