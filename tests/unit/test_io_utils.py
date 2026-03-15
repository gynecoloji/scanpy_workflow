import pytest
import anndata as ad
import numpy as np
import scipy.sparse
from scanpy_workflow.utils.io import read_h5ad, write_h5ad


def test_write_and_read_roundtrip(tmp_path, small_adata):
    out = tmp_path / "test.h5ad"
    write_h5ad(small_adata, str(out))
    loaded = read_h5ad(str(out))
    assert loaded.shape == small_adata.shape
    assert list(loaded.obs_names) == list(small_adata.obs_names)


def test_write_creates_parent_dirs(tmp_path, small_adata):
    out = tmp_path / "subdir" / "nested" / "test.h5ad"
    write_h5ad(small_adata, str(out))
    assert out.exists()
