import pytest
from scanpy_workflow.preprocessing.hvg import select_hvg


def test_hvg_seurat_v3_uses_counts_layer(normalized_adata):
    adata = select_hvg(normalized_adata, n_top_genes=50, flavor="seurat_v3")
    assert "highly_variable" in adata.var.columns
    assert adata.var["highly_variable"].sum() <= 50


def test_hvg_seurat_v3_requires_counts_layer(small_adata):
    with pytest.raises(ValueError, match="layers\\['counts'\\] required"):
        select_hvg(small_adata, n_top_genes=50, flavor="seurat_v3")


def test_hvg_seurat_flavor_uses_x(normalized_adata):
    adata = select_hvg(normalized_adata, n_top_genes=50, flavor="seurat")
    assert "highly_variable" in adata.var.columns


def test_hvg_marks_correct_count(normalized_adata):
    n = 30
    adata = select_hvg(normalized_adata.copy(), n_top_genes=n, flavor="seurat")
    assert adata.var["highly_variable"].sum() <= n
