"""
CLI entrypoints for scanpy_workflow.
All multi-word subcommands use hyphens (Unix convention).
R-backed steps (ambient, normalize scran, impute alra, annotate sctype)
are called directly by Nextflow via Rscript — they are NOT registered here.
"""
import click
from scanpy_workflow.utils.io import read_h5ad, write_h5ad


@click.group()
def cli():
    """Scanpy end-to-end single-cell workflow."""
    pass


# ── load ──────────────────────────────────────────────────────────────────────
@cli.command("load")
@click.option("--input", "input_dir", required=True, help="Cell Ranger output directory")
@click.option("--output", "output_path", required=True, help="Output .h5ad path")
def cmd_load(input_dir, output_path):
    """Load Cell Ranger output (MEX or HDF5)."""
    from scanpy_workflow.io.load import load_cellranger
    adata = load_cellranger(input_dir)
    write_h5ad(adata, output_path)
    click.echo(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes -> {output_path}")


# ── qc ────────────────────────────────────────────────────────────────────────
@cli.command("qc")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--plots-dir", "plots_dir", default=None)
def cmd_qc(input_path, output_path, plots_dir):
    """Compute QC metrics."""
    from scanpy_workflow.qc.metrics import compute_qc_metrics
    from scanpy_workflow.qc.plots import plot_qc
    adata = read_h5ad(input_path)
    adata = compute_qc_metrics(adata)
    if plots_dir:
        plot_qc(adata, output_dir=plots_dir)
    write_h5ad(adata, output_path)


# ── qc-report ─────────────────────────────────────────────────────────────────
@cli.command("qc-report")
@click.option("--input",        "input_path",  required=True)
@click.option("--sample",       "sample_name", default="sample",
              help="Sample label for the report header [default: sample]")
@click.option("--output",       "output_path", required=True,
              help="Output .html path")
@click.option("--min-genes",    default=200,   type=int)
@click.option("--max-genes",    default=6000,  type=int)
@click.option("--min-counts",   default=500,   type=int)
@click.option("--max-counts",   default=30000, type=int)
@click.option("--max-pct-mito", default=20.0,  type=float)
@click.option("--max-pct-ribo", default=50.0,  type=float)
def cmd_qc_report(input_path, sample_name, output_path,
                  min_genes, max_genes, min_counts, max_counts,
                  max_pct_mito, max_pct_ribo):
    """Generate per-sample QC HTML report with filter-threshold overlays."""
    from scanpy_workflow.qc.metrics import compute_qc_metrics
    from scanpy_workflow.qc.report import generate_qc_report
    adata = read_h5ad(input_path)
    if "n_genes_by_counts" not in adata.obs.columns:
        adata = compute_qc_metrics(adata)
    thresholds = dict(
        min_genes=min_genes, max_genes=max_genes,
        min_counts=min_counts, max_counts=max_counts,
        max_pct_mito=max_pct_mito, max_pct_ribo=max_pct_ribo,
    )
    generate_qc_report(adata, sample_name=sample_name,
                       thresholds=thresholds, output_path=output_path)
    click.echo(f"QC report written to {output_path}")


# ── filter ────────────────────────────────────────────────────────────────────
@cli.command("filter")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--min-genes", default=200, type=int)
@click.option("--max-genes", default=6000, type=int)
@click.option("--min-counts", default=500, type=int)
@click.option("--max-counts", default=30000, type=int)
@click.option("--max-pct-mito", default=20.0, type=float)
@click.option("--max-pct-ribo", default=50.0, type=float)
@click.option("--min-cells", default=3, type=int)
def cmd_filter(input_path, output_path, min_genes, max_genes, min_counts,
               max_counts, max_pct_mito, max_pct_ribo, min_cells):
    """Filter cells and genes by QC cutoffs."""
    from scanpy_workflow.qc.filter import filter_cells_genes
    adata = read_h5ad(input_path)
    adata = filter_cells_genes(adata, min_genes=min_genes, max_genes=max_genes,
                                min_counts=min_counts, max_counts=max_counts,
                                max_pct_mito=max_pct_mito, max_pct_ribo=max_pct_ribo,
                                min_cells=min_cells)
    write_h5ad(adata, output_path)
    click.echo(f"Filtered: {adata.n_obs} cells, {adata.n_vars} genes -> {output_path}")


# ── doublets ──────────────────────────────────────────────────────────────────
@cli.command("doublets")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--filter-doublets", is_flag=True, default=False)
def cmd_doublets(input_path, output_path, filter_doublets):
    """Detect (and optionally filter) doublets with Scrublet."""
    from scanpy_workflow.doublets.scrublet import detect_doublets
    adata = read_h5ad(input_path)
    adata = detect_doublets(adata, filter_doublets=filter_doublets)
    write_h5ad(adata, output_path)


# ── normalize ─────────────────────────────────────────────────────────────────
@cli.command("normalize")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", default="library_size", type=click.Choice(["library_size"]))
@click.option("--target-sum", default=10000.0, type=float)
def cmd_normalize(input_path, output_path, method, target_sum):
    """Normalize counts (library_size). For scran, use scran.R."""
    from scanpy_workflow.preprocessing.normalize import normalize
    adata = read_h5ad(input_path)
    adata = normalize(adata, method=method, target_sum=target_sum)
    write_h5ad(adata, output_path)


# ── hvg ───────────────────────────────────────────────────────────────────────
@cli.command("hvg")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--n-top-genes", default=3000, type=int)
@click.option("--flavor", default="seurat_v3",
              type=click.Choice(["seurat_v3", "seurat", "cell_ranger"]))
def cmd_hvg(input_path, output_path, n_top_genes, flavor):
    """Select highly variable genes."""
    from scanpy_workflow.preprocessing.hvg import select_hvg
    adata = read_h5ad(input_path)
    adata = select_hvg(adata, n_top_genes=n_top_genes, flavor=flavor)
    write_h5ad(adata, output_path)


# ── merge ─────────────────────────────────────────────────────────────────────
@cli.command("merge")
@click.option("--inputs", "input_paths", required=True, multiple=True,
              help="Paths to per-sample .h5ad files")
@click.option("--sample-names", required=True, multiple=True)
@click.option("--output", "output_path", required=True)
@click.option("--batch-key", default="sample")
def cmd_merge(input_paths, sample_names, output_path, batch_key):
    """Merge per-sample AnnData objects."""
    from scanpy_workflow.preprocessing.merge import merge_samples
    adatas = [read_h5ad(p) for p in input_paths]
    merged = merge_samples(adatas, sample_names=list(sample_names), batch_key=batch_key)
    write_h5ad(merged, output_path)
    click.echo(f"Merged {len(adatas)} samples: {merged.n_obs} cells -> {output_path}")


# ── scale ─────────────────────────────────────────────────────────────────────
@cli.command("scale")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--max-value", default=10.0, type=float)
def cmd_scale(input_path, output_path, max_value):
    """Z-score scale gene expression."""
    from scanpy_workflow.preprocessing.scale import scale
    adata = read_h5ad(input_path)
    adata = scale(adata, max_value=max_value)
    write_h5ad(adata, output_path)


# ── pca ───────────────────────────────────────────────────────────────────────
@cli.command("pca")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--n-comps", default=50, type=int)
@click.option("--scvi-path", is_flag=True, default=False,
              help="Assert X equals layers['log_norm'] (scVI path, scaling was skipped)")
def cmd_pca(input_path, output_path, n_comps, scvi_path):
    """Run PCA."""
    from scanpy_workflow.preprocessing.pca import run_pca
    adata = read_h5ad(input_path)
    adata = run_pca(adata, n_comps=n_comps, scvi_path=scvi_path)
    write_h5ad(adata, output_path)


# ── batch-correct ──────────────────────────────────────────────────────────────
@cli.command("batch-correct")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", required=True, type=click.Choice(["harmony", "bbknn", "scvi"]))
@click.option("--batch-key", default="sample")
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
def cmd_batch_correct(input_path, output_path, method, batch_key, evaluate, eval_dir):
    """Apply batch correction."""
    adata = read_h5ad(input_path)
    if method == "harmony":
        from scanpy_workflow.integration.harmony import batch_correct_harmony
        adata = batch_correct_harmony(adata, batch_key=batch_key,
                                       run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "bbknn":
        from scanpy_workflow.integration.bbknn import batch_correct_bbknn
        adata = batch_correct_bbknn(adata, batch_key=batch_key,
                                     run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "scvi":
        from scanpy_workflow.integration.scvi import batch_correct_scvi
        adata = batch_correct_scvi(adata, batch_key=batch_key,
                                    run_evaluation=evaluate, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── impute ────────────────────────────────────────────────────────────────────
@cli.command("impute")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--method", required=True, type=click.Choice(["magic", "scvi"]))
@click.option("--batch-key", default="sample")
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
@click.option("--model-dir", default=None)
def cmd_impute(input_path, output_path, method, batch_key, evaluate, eval_dir, model_dir):
    """Impute missing values (Python methods: magic, scvi). For alra, use alra.R."""
    adata = read_h5ad(input_path)
    if method == "magic":
        from scanpy_workflow.imputation.magic import impute_magic
        adata = impute_magic(adata, run_evaluation=evaluate, output_dir=eval_dir)
    elif method == "scvi":
        from scanpy_workflow.imputation.scvi import impute_scvi
        adata = impute_scvi(adata, batch_key=batch_key, model_dir=model_dir,
                             run_evaluation=evaluate, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── cluster ───────────────────────────────────────────────────────────────────
@cli.command("cluster")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--algorithm", default="leiden", type=click.Choice(["leiden", "louvain"]))
@click.option("--resolution", default=0.5, type=float)
@click.option("--n-neighbors", default=15, type=int)
@click.option("--evaluate", is_flag=True, default=False)
@click.option("--eval-dir", default=None)
def cmd_cluster(input_path, output_path, algorithm, resolution, n_neighbors, evaluate, eval_dir):
    """Build neighbor graph and cluster cells."""
    from scanpy_workflow.clustering.cluster import cluster
    adata = read_h5ad(input_path)
    adata = cluster(adata, algorithm=algorithm, resolution=resolution, n_neighbors=n_neighbors)
    if evaluate and eval_dir:
        from scanpy_workflow.clustering.evaluate import evaluate_clustering
        evaluate_clustering(adata, cluster_key=algorithm, output_dir=eval_dir)
    write_h5ad(adata, output_path)


# ── rank-genes ────────────────────────────────────────────────────────────────
@cli.command("rank-genes")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--groupby", default="leiden")
@click.option("--output-dir", "marker_dir", default=None)
def cmd_rank_genes(input_path, output_path, groupby, marker_dir):
    """Identify marker genes per cluster."""
    from scanpy_workflow.annotation.markers import rank_marker_genes
    adata = read_h5ad(input_path)
    adata = rank_marker_genes(adata, groupby=groupby, output_dir=marker_dir)
    write_h5ad(adata, output_path)


# ── annotate ──────────────────────────────────────────────────────────────────
@cli.command("annotate")
@click.option("--input", "input_path", required=True)
@click.option("--output", "output_path", required=True)
@click.option("--model", default="Immune_All_Low.pkl")
@click.option("--output-dir", "annot_dir", default=None)
def cmd_annotate(input_path, output_path, model, annot_dir):
    """Automated annotation with CellTypist."""
    from scanpy_workflow.annotation.celltypist import annotate_celltypist
    adata = read_h5ad(input_path)
    adata = annotate_celltypist(adata, model=model, output_dir=annot_dir)
    write_h5ad(adata, output_path)


# ── de ────────────────────────────────────────────────────────────────────────
@cli.command("de")
@click.option("--input", "input_path", required=True)
@click.option("--output-h5ad", "output_path", required=True)
@click.option("--output-dir", "de_dir", required=True)
@click.option("--groupby", default="leiden")
@click.option("--method", default="wilcoxon",
              type=click.Choice(["wilcoxon", "t-test", "logreg"]))
def cmd_de(input_path, output_path, de_dir, groupby, method):
    """Run differential expression."""
    from scanpy_workflow.de.de import run_de
    adata = read_h5ad(input_path)
    run_de(adata, groupby=groupby, method=method, output_dir=de_dir)
    write_h5ad(adata, output_path)


# ── pathway-score ─────────────────────────────────────────────────────────────
@cli.command("pathway-score")
@click.option("--input", "input_path", required=True,
              help="Input .h5ad (log-normalised X)")
@click.option("--output", "output_path", required=True,
              help="Output .h5ad with pathway scores in obsm['pathway_scores']")
@click.option("--source", required=True,
              type=click.Choice([
                  "progeny", "msigdb_hallmark", "msigdb_kegg",
                  "msigdb_reactome", "msigdb_gobp", "custom"
              ]),
              help="Gene-set source")
@click.option("--method", default="aucell",
              type=click.Choice(["aucell", "ulm", "scanpy_score"]),
              help="Scoring method [default: aucell]")
@click.option("--groupby", default="leiden",
              help="obs column for per-cluster summary [default: leiden]")
@click.option("--min-n", "min_n", default=5, type=int,
              help="Min gene-set / data overlap to keep a gene set [default: 5]")
@click.option("--organism", default="human",
              type=click.Choice(["human", "mouse"]),
              help="Organism for built-in gene sets [default: human]")
@click.option("--custom-genesets", "custom_path", default=None,
              help="Path to .gmt or .csv gene-set file (required when source=custom)")
@click.option("--output-dir", "output_dir", default=None,
              help="Directory to write pathway_scores.csv and cluster summary")
def cmd_pathway_score(input_path, output_path, source, method, groupby,
                      min_n, organism, custom_path, output_dir):
    """Score cells against pathway / gene-set collections."""
    from scanpy_workflow.pathway.genesets import load_genesets
    from scanpy_workflow.pathway.score import run_pathway_scoring
    adata = read_h5ad(input_path)
    net = load_genesets(source, organism=organism, custom_path=custom_path)
    scores = run_pathway_scoring(
        adata, net, method=method, groupby=groupby,
        min_n=min_n, output_dir=output_dir
    )
    write_h5ad(adata, output_path)
    click.echo(
        f"Scored {scores.shape[1]} gene sets across {scores.shape[0]} cells -> {output_path}"
    )


# ── report ────────────────────────────────────────────────────────────────────
@cli.command("report")
@click.option("--input", "input_path", required=True)
@click.option("--output-dir", "output_dir", required=True)
@click.option("--config", "config_path", default=None)
@click.option("--cluster-key", default="leiden")
def cmd_report(input_path, output_dir, config_path, cluster_key):
    """Generate HTML summary report."""
    from scanpy_workflow.reporting.report import generate_report
    from scanpy_workflow.utils.config import load_config
    adata = read_h5ad(input_path)
    params = load_config(config_path) if config_path else {}
    generate_report(adata, output_dir=output_dir, params=params, cluster_key=cluster_key)
    click.echo(f"Report written to {output_dir}/report.html")
