// workflow/modules/doublets.nf
process DOUBLETS {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/doublets", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_doublets.h5ad")

    script:
    """
    scanpy-workflow doublets \
        --input  ${h5ad} \
        --output ${sample}_doublets.h5ad
    """
}
