import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import anndata as ad
from scanpy_workflow.pathway.genesets import load_genesets, _parse_gmt, _load_custom_file
from scanpy_workflow.pathway.score import run_pathway_scoring

REQUIRED_COLS = ["source", "target"]


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def small_adata():
    """200-cell, 80-gene AnnData with log-normalised X and leiden labels."""
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(200, 80)).astype(float)
    # log-normalise in place
    totals = X.sum(axis=1, keepdims=True)
    X = np.log1p(X / totals * 1e4)
    gene_names = [f"Gene{i}" for i in range(80)]
    cell_names = [f"Cell{i}" for i in range(200)]
    adata = ad.AnnData(X, obs=pd.DataFrame(index=cell_names),
                       var=pd.DataFrame(index=gene_names))
    adata.obs["leiden"] = np.tile(["c1", "c2", "c3", "c4"], 50)
    return adata


@pytest.fixture
def small_net():
    """Minimal gene-set network: 3 pathways × 15 genes each (from Gene0–Gene44)."""
    rows = []
    for i, gs in enumerate(["PathwayA", "PathwayB", "PathwayC"]):
        for g in [f"Gene{j}" for j in range(i * 15, (i + 1) * 15)]:
            rows.append({"source": gs, "target": g, "weight": 1.0})
    return pd.DataFrame(rows)


@pytest.fixture
def gmt_file(tmp_path):
    gmt = tmp_path / "test.gmt"
    gmt.write_text(
        "PathwayA\tA desc\tGene0\tGene1\tGene2\tGene3\tGene4\n"
        "PathwayB\tB desc\tGene5\tGene6\tGene7\tGene8\tGene9\n"
    )
    return str(gmt)


@pytest.fixture
def csv_file(tmp_path):
    rows = [{"gene_set": f"Path{g}", "gene": f"Gene{i}"}
            for g in ["X", "Y"] for i in range(10)]
    csv = tmp_path / "custom.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    return str(csv)


# ── genesets tests ─────────────────────────────────────────────────────────────

def test_parse_gmt_returns_source_target(gmt_file):
    df = load_genesets("custom", custom_path=gmt_file)
    assert set(df.columns) >= {"source", "target"}
    assert set(df["source"].unique()) == {"PathwayA", "PathwayB"}
    assert df["target"].str.startswith("Gene").all()


def test_parse_gmt_gene_count(gmt_file):
    df = load_genesets("custom", custom_path=gmt_file)
    assert len(df) == 10  # 5 + 5 genes


def test_csv_custom_genesets(csv_file):
    df = load_genesets("custom", custom_path=csv_file)
    assert set(df.columns) >= {"source", "target"}
    assert set(df["source"].unique()) == {"PathX", "PathY"}


def test_unknown_source_raises():
    with pytest.raises(ValueError, match="Unknown gene-set source"):
        load_genesets("not_a_source")


def test_custom_without_path_raises():
    with pytest.raises(ValueError, match="custom_path is required"):
        load_genesets("custom")


def test_missing_custom_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_genesets("custom", custom_path=str(tmp_path / "nonexistent.gmt"))


# ── scoring tests ─────────────────────────────────────────────────────────────

def test_scanpy_score_returns_dataframe(small_adata, small_net):
    scores = run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden", min_n=5
    )
    assert isinstance(scores, pd.DataFrame)
    assert scores.shape == (200, 3)  # 200 cells × 3 pathways


def test_scanpy_score_columns_match_pathways(small_adata, small_net):
    scores = run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden", min_n=5
    )
    assert set(scores.columns) == {"PathwayA", "PathwayB", "PathwayC"}


def test_scanpy_score_stored_in_obsm(small_adata, small_net):
    run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden", min_n=5
    )
    assert "pathway_scores" in small_adata.obsm
    assert isinstance(small_adata.obsm["pathway_scores"], pd.DataFrame)


def test_scanpy_score_no_obs_column_pollution(small_adata, small_net):
    obs_cols_before = set(small_adata.obs.columns)
    run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden", min_n=5
    )
    assert set(small_adata.obs.columns) == obs_cols_before


def test_scanpy_score_writes_csvs(tmp_path, small_adata, small_net):
    run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden",
        min_n=5, output_dir=str(tmp_path)
    )
    assert (tmp_path / "pathway_scores.csv").exists()
    assert (tmp_path / "cluster_pathway_scores.csv").exists()


def test_cluster_summary_shape(tmp_path, small_adata, small_net):
    run_pathway_scoring(
        small_adata, small_net, method="scanpy_score", groupby="leiden",
        min_n=5, output_dir=str(tmp_path)
    )
    summary = pd.read_csv(tmp_path / "cluster_pathway_scores.csv", index_col=0)
    assert summary.shape == (4, 3)  # 4 clusters × 3 pathways


def test_min_n_filters_small_genesets(small_adata):
    """Gene sets with fewer than min_n matching genes should be dropped."""
    tiny_net = pd.DataFrame([
        {"source": "TinyPath", "target": "Gene0"},
        {"source": "TinyPath", "target": "Gene1"},
        {"source": "BigPath",  "target": f"Gene{i}"} for i in range(20)
    ])
    scores = run_pathway_scoring(
        small_adata, tiny_net, method="scanpy_score", min_n=5
    )
    assert "TinyPath" not in scores.columns
    assert "BigPath" in scores.columns


def test_invalid_method_raises(small_adata, small_net):
    with pytest.raises(ValueError, match="Unknown method"):
        run_pathway_scoring(small_adata, small_net, method="bad_method")


def test_missing_net_columns_raises(small_adata):
    bad_net = pd.DataFrame([{"pathway": "A", "gene": "Gene0"}])
    with pytest.raises(ValueError, match="source.*target"):
        run_pathway_scoring(small_adata, bad_net, method="scanpy_score")


@pytest.mark.skipif(
    pytest.importorskip("decoupler", reason="decoupler not installed") is None,
    reason="decoupler not installed"
)
def test_aucell_returns_dataframe(small_adata, small_net):
    pytest.importorskip("decoupler")
    scores = run_pathway_scoring(
        small_adata, small_net, method="aucell", groupby="leiden", min_n=5
    )
    assert isinstance(scores, pd.DataFrame)
    assert scores.shape[0] == 200


@pytest.mark.skipif(
    pytest.importorskip("decoupler", reason="decoupler not installed") is None,
    reason="decoupler not installed"
)
def test_aucell_obsm_cleanup(small_adata, small_net):
    """Intermediate obsm keys (aucell_estimate, aucell_pvals) must be removed."""
    pytest.importorskip("decoupler")
    run_pathway_scoring(
        small_adata, small_net, method="aucell", groupby="leiden", min_n=5
    )
    assert "aucell_estimate" not in small_adata.obsm
    assert "aucell_pvals" not in small_adata.obsm
