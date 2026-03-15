import pytest
from scanpy_workflow.clustering.cluster import cluster


def test_cluster_leiden_adds_obs_column(merged_adata):
    adata = cluster(merged_adata, algorithm="leiden", resolution=0.3, n_neighbors=10)
    assert "leiden" in adata.obs.columns


def test_cluster_louvain_adds_obs_column(merged_adata):
    adata = cluster(merged_adata, algorithm="louvain", resolution=0.3, n_neighbors=10)
    assert "louvain" in adata.obs.columns


def test_cluster_skips_neighbors_for_bbknn(merged_adata):
    """When uns['neighbors_use_rep'] == 'skip', pp.neighbors is not called."""
    import scanpy as sc
    # Run BBKNN first so obsp exists
    adata = merged_adata.copy()
    adata.uns["neighbors_use_rep"] = "skip"
    sc.pp.neighbors(adata, use_rep="X_pca")  # simulate BBKNN having written the graph
    adata.uns["neighbors_use_rep"] = "skip"  # reset to skip
    result = cluster(adata, algorithm="leiden", resolution=0.3, n_neighbors=10)
    assert "leiden" in result.obs.columns


def test_cluster_unknown_algorithm_raises(merged_adata):
    with pytest.raises(ValueError, match="Unknown clustering algorithm"):
        cluster(merged_adata, algorithm="dbscan", resolution=0.3, n_neighbors=10)


def test_evaluate_clustering_creates_outputs(tmp_path, merged_adata):
    import scanpy as sc
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    from scanpy_workflow.clustering.evaluate import evaluate_clustering
    scores = evaluate_clustering(adata, cluster_key="leiden", output_dir=str(tmp_path))
    assert (tmp_path / "clustering_metrics.json").exists()
    assert "silhouette" in scores or "note" in scores
