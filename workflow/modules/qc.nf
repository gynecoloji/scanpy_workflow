// workflow/modules/qc.nf
process QC {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/qc", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_qc.h5ad")
    path  "qc_plots/"

    script:
    """
    scanpy-workflow qc \
        --input  ${h5ad} \
        --output ${sample}_qc.h5ad \
        --plots  qc_plots/
    """
}
