import matplotlib
matplotlib.use("Agg")
from pathlib import Path
from scanpy_workflow.qc.metrics import compute_qc_metrics
from scanpy_workflow.qc.plots import plot_qc


def test_plot_qc_creates_files(tmp_path, small_adata):
    adata = compute_qc_metrics(small_adata)
    plot_qc(adata, output_dir=str(tmp_path))
    assert (tmp_path / "qc_violin.png").exists()
    assert (tmp_path / "qc_scatter_counts_vs_genes.png").exists()
