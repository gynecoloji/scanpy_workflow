// workflow/modules/annotation.nf
process RANK_GENES {
    label     'scanpy'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path h5ad
    val  groupby

    output:
    path "rank_genes.h5ad"
    path "marker_plots/"

    script:
    """
    scanpy-workflow rank-genes \
        --input   ${h5ad} \
        --output  rank_genes.h5ad \
        --groupby ${groupby} \
        --plots   marker_plots/
    """
}

process ANNOTATE_CELLTYPIST {
    label     'scanpy'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path  h5ad
    val   model
    val   groupby

    output:
    path "annotated_celltypist.h5ad"

    script:
    """
    scanpy-workflow annotate \
        --input   ${h5ad} \
        --output  annotated_celltypist.h5ad \
        --model   ${model} \
        --groupby ${groupby}
    """
}

process ANNOTATE_SCTYPE {
    label     'r_env'
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    path h5ad
    path markers_json
    val  groupby

    output:
    path "annotated_sctype.h5ad"

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/annotation/sctype.R \
        --input   ${h5ad} \
        --markers ${markers_json} \
        --groupby ${groupby} \
        --output  annotated_sctype.h5ad
    """
}
