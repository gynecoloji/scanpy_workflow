// workflow/modules/load.nf
process LOAD {
    label 'scanpy'
    tag   "${sample}"

    input:
    tuple val(sample), path(input_dir)

    output:
    tuple val(sample), path("${sample}_loaded.h5ad")

    script:
    """
    scanpy-workflow load \
        --input    ${input_dir} \
        --sample   ${sample} \
        --output   ${sample}_loaded.h5ad
    """
}
