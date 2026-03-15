// workflow/modules/scale.nf
process SCALE {
    label     'scanpy'
    publishDir "${params.output_dir}/scaled", mode: 'copy'

    input:
    path h5ad

    output:
    path "scaled.h5ad"

    script:
    """
    scanpy-workflow scale \
        --input  ${h5ad} \
        --output scaled.h5ad
    """
}
