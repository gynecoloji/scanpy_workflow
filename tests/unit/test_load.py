import pytest
import numpy as np
from pathlib import Path
from scanpy_workflow.io.load import load_cellranger


def test_load_cellranger_mex(tmp_path):
    """Test loading MEX format (barcodes/features/matrix)."""
    import scanpy as sc
    import scipy.io
    import gzip
    import anndata as ad

    # Create minimal MEX structure
    sample_dir = tmp_path / "sample_mex" / "filtered_feature_bc_matrix"
    sample_dir.mkdir(parents=True)

    n_cells, n_genes = 50, 100
    rng = np.random.default_rng(0)
    matrix = scipy.sparse.random(n_genes, n_cells, density=0.1, format="csc",
                                  random_state=0, data_rvs=lambda s: rng.integers(1, 10, s))
    # Write matrix as gzipped .mtx to match scanpy's expected format
    import io, gzip as gz
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, matrix)
    with gz.open(sample_dir / "matrix.mtx.gz", "wb") as f:
        f.write(buf.getvalue())

    with gzip.open(sample_dir / "barcodes.tsv.gz", "wt") as f:
        for i in range(n_cells):
            f.write(f"CELL{i:04d}-1\n")

    with gzip.open(sample_dir / "features.tsv.gz", "wt") as f:
        for i in range(n_genes):
            f.write(f"ENSG{i:08d}\tGene{i}\tGene Expression\n")

    adata = load_cellranger(str(tmp_path / "sample_mex"))
    assert adata.shape == (n_cells, n_genes)
    assert adata.X is not None


def test_load_raises_if_no_data(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        load_cellranger(str(empty_dir))
