// workflow/modules/report.nf
process REPORT {
    label     'scanpy'
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
    path h5ad
    path de_csv
    val  output_dir

    output:
    path "report.html"
    path "final.h5ad"

    script:
    """
    scanpy-workflow report \
        --input      ${h5ad} \
        --de-csv     ${de_csv} \
        --output-dir . \
        --output     report.html \
        --final-h5ad final.h5ad
    """
}
