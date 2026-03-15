import scanpy as sc
from scanpy_workflow.annotation.markers import rank_marker_genes


def test_rank_genes_adds_uns_key(merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    result = rank_marker_genes(adata, groupby="leiden")
    assert "rank_genes_groups" in result.uns


def test_rank_genes_csv_written(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    rank_marker_genes(adata, groupby="leiden", output_dir=str(tmp_path))
    assert (tmp_path / "marker_genes.csv").exists()
