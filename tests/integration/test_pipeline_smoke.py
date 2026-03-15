# tests/integration/test_pipeline_smoke.py
"""
Smoke test: run the Nextflow pipeline on a synthetic 2-sample dataset.
Verifies:
  - Pipeline exits 0
  - Key output files exist (final.h5ad, de_results.csv, report.html)
  - final.h5ad has the expected AnnData structure
Skipped if Nextflow is not on PATH (CI installs it; local dev may skip).
"""
import os
import subprocess
import shutil
import pytest
import anndata as ad


def nextflow_available():
    return shutil.which("nextflow") is not None


@pytest.mark.skipif(not nextflow_available(), reason="Nextflow not on PATH")
def test_pipeline_smoke(cellranger_dir, tmp_path):
    out_dir = str(tmp_path / "results")
    params  = str(tmp_path / "params.yaml")

    # Write minimal params.yaml for smoke run
    with open(params, "w") as f:
        f.write(f"""\
input_dir: "{cellranger_dir}"
samples: [sample_A, sample_B]
output_dir: "{out_dir}"
ambient:
  enabled: false
doublet:
  enabled: false
normalization:
  method: library_size
  target_sum: 10000
hvg:
  n_top_genes: 50
  post_merge_n_top_genes: 30
  flavor: seurat_v3
pca:
  n_comps: 10
batch_correction:
  enabled: true
  method: harmony
  batch_key: sample
imputation:
  enabled: false
clustering:
  algorithm: leiden
  resolution: 0.5
  n_neighbors: 5
  evaluate: false
annotation:
  unsupervised: true
  automated: false
de:
  method: wilcoxon
  groupby: leiden
report: true
""")

    # Resolve repo root so relative paths in main.nf work regardless of pytest invocation dir
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    result = subprocess.run(
        ["nextflow", "run", "workflow/main.nf",
         "-params-file", params,
         "-profile", "conda",
         "-resume"],
        capture_output=True, text=True, timeout=600,
        cwd=repo_root
    )

    assert result.returncode == 0, (
        f"Nextflow pipeline failed:\nSTDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Key output files
    assert os.path.exists(os.path.join(out_dir, "report", "final.h5ad")), \
        "Missing final.h5ad"
    assert os.path.exists(os.path.join(out_dir, "de", "de_results.csv")), \
        "Missing de_results.csv"
    assert os.path.exists(os.path.join(out_dir, "report", "report.html")), \
        "Missing report.html"

    # Structural checks on final.h5ad
    adata = ad.read_h5ad(os.path.join(out_dir, "report", "final.h5ad"))
    assert "leiden" in adata.obs.columns, "Missing leiden column in obs"
    assert "counts" in adata.layers,      "Missing counts layer"
    assert "log_norm" in adata.layers,    "Missing log_norm layer"

    # DE CSV has required columns
    import pandas as pd
    de = pd.read_csv(os.path.join(out_dir, "de", "de_results.csv"))
    required_cols = {"group", "gene", "score", "logfoldchange", "pval", "pval_adj"}
    assert required_cols.issubset(set(de.columns)), \
        f"DE CSV missing columns: {required_cols - set(de.columns)}"
