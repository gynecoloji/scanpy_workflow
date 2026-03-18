// workflow/modules/qc_report.nf
process QC_REPORT {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/qc_report", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    path "qc_report.html"

    script:
    """
    scanpy-workflow qc-report \
        --input        ${h5ad} \
        --sample       "${sample}" \
        --output       qc_report.html \
        --min-genes    ${params.filter.min_genes} \
        --max-genes    ${params.filter.max_genes} \
        --min-counts   ${params.filter.min_counts} \
        --max-counts   ${params.filter.max_counts} \
        --max-pct-mito ${params.filter.max_pct_mito} \
        --max-pct-ribo ${params.filter.max_pct_ribo}
    """
}
