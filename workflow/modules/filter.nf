// workflow/modules/filter.nf
process FILTER {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/filtered", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_filtered.h5ad")

    script:
    """
    scanpy-workflow filter \
        --input        ${h5ad} \
        --output       ${sample}_filtered.h5ad \
        --min-genes    ${params.filter.min_genes} \
        --max-genes    ${params.filter.max_genes} \
        --min-counts   ${params.filter.min_counts} \
        --max-counts   ${params.filter.max_counts} \
        --max-pct-mito ${params.filter.max_pct_mito} \
        --max-pct-ribo ${params.filter.max_pct_ribo} \
        --min-cells    ${params.filter.min_cells}
    """
}
