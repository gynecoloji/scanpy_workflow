"""
End-to-end smoke test using a programmatically generated synthetic dataset.
No external data downloads required. Tests the full Python package pipeline.
"""
import pytest
import numpy as np
import anndata as ad
import scipy.sparse
import scanpy as sc
from pathlib import Path


@pytest.fixture(scope="module")
def synthetic_cellranger(tmp_path_factory):
    """Create two synthetic Cell Ranger MEX directories."""
    import scipy.io, gzip
    base = tmp_path_factory.mktemp("cellranger")
    sample_dirs = []
    for s in range(2):
        rng = np.random.default_rng(s)
        n_cells, n_genes = 100, 200
        sample_dir = base / f"sample_{s}" / "filtered_feature_bc_matrix"
        sample_dir.mkdir(parents=True)
        matrix = scipy.sparse.random(n_genes, n_cells, density=0.3, format="csc",
                                      random_state=s, data_rvs=lambda n: rng.integers(1, 20, n))
        import io, gzip as gz
        buf = io.BytesIO()
        scipy.io.mmwrite(buf, matrix)
        with gz.open(sample_dir / "matrix.mtx.gz", "wb") as f:
            f.write(buf.getvalue())
        with gzip.open(sample_dir / "barcodes.tsv.gz", "wt") as f:
            for i in range(n_cells):
                f.write(f"CELL{s}{i:04d}-1\n")
        with gzip.open(sample_dir / "features.tsv.gz", "wt") as f:
            for i in range(n_genes):
                prefix = "MT-" if i < 5 else "RPS" if i < 10 else "Gene"
                f.write(f"ENSG{i:08d}\t{prefix}{i}\tGene Expression\n")
        sample_dirs.append(str(sample_dir.parent))
    return sample_dirs


def test_full_pipeline_smoke(synthetic_cellranger, tmp_path):
    from scanpy_workflow.io.load import load_cellranger
    from scanpy_workflow.qc.metrics import compute_qc_metrics
    from scanpy_workflow.qc.filter import filter_cells_genes
    from scanpy_workflow.doublets.scrublet import detect_doublets
    from scanpy_workflow.preprocessing.normalize import normalize
    from scanpy_workflow.preprocessing.hvg import select_hvg
    from scanpy_workflow.preprocessing.merge import merge_samples
    from scanpy_workflow.preprocessing.scale import scale
    from scanpy_workflow.preprocessing.pca import run_pca
    from scanpy_workflow.integration.harmony import batch_correct_harmony
    from scanpy_workflow.clustering.cluster import cluster
    from scanpy_workflow.de.de import run_de

    # Per-sample
    adatas = []
    for s, sample_dir in enumerate(synthetic_cellranger):
        adata = load_cellranger(sample_dir)
        adata = compute_qc_metrics(adata)
        adata = filter_cells_genes(adata, min_genes=1, max_genes=10000,
                                    min_counts=1, max_counts=10**9,
                                    max_pct_mito=100, max_pct_ribo=100, min_cells=1)
        adata = detect_doublets(adata, filter_doublets=False)
        adata = normalize(adata)
        adata = select_hvg(adata, n_top_genes=50, flavor="seurat_v3")
        adatas.append(adata)

    # Merge + shared
    merged = merge_samples(adatas, sample_names=["s0", "s1"], batch_key="sample")
    select_hvg(merged, n_top_genes=50, flavor="seurat_v3")
    scale(merged)
    run_pca(merged, n_comps=10)
    batch_correct_harmony(merged, batch_key="sample")
    cluster(merged, algorithm="leiden", resolution=0.5, n_neighbors=5)
    df = run_de(merged, groupby="leiden", method="wilcoxon")

    # Assertions
    assert merged.n_obs > 0
    assert "leiden" in merged.obs.columns
    assert "X_pca_harmony" in merged.obsm
    assert set(df.columns) == {"group", "gene", "score", "logfoldchange", "pval", "pval_adj"}
