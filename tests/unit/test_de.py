import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scanpy_workflow.de.de import run_de

REQUIRED_COLS = ["group", "gene", "score", "logfoldchange", "pval", "pval_adj"]


def _clustered_adata(merged_adata):
    adata = merged_adata.copy()
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=0.3)
    return adata


def test_de_returns_required_columns(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="wilcoxon")
    for col in REQUIRED_COLS:
        assert col in df.columns, f"Missing column: {col}"


def test_de_logreg_pval_is_nan(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="logreg")
    assert df["pval"].isna().all()
    assert df["pval_adj"].isna().all()


def test_de_wilcoxon_pval_not_nan(merged_adata):
    adata = _clustered_adata(merged_adata)
    df = run_de(adata, groupby="leiden", method="wilcoxon")
    assert not df["pval"].isna().all()


def test_de_raises_for_missing_groupby(merged_adata):
    with pytest.raises(ValueError, match="groupby 'nonexistent' not found"):
        run_de(merged_adata, groupby="nonexistent", method="wilcoxon")


def test_de_csv_written(tmp_path, merged_adata):
    adata = _clustered_adata(merged_adata)
    run_de(adata, groupby="leiden", method="wilcoxon", output_dir=str(tmp_path))
    assert (tmp_path / "de_results.csv").exists()
