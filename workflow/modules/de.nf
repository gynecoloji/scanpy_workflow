// workflow/modules/de.nf
process DE {
    label     'scanpy'
    publishDir "${params.output_dir}/de", mode: 'copy'

    input:
    path h5ad
    val  method    // "wilcoxon" | "t-test" | "logreg"
    val  groupby

    output:
    path "de_results.csv"
    path "de_plots/"

    script:
    """
    scanpy-workflow de \
        --input   ${h5ad} \
        --output  de_results.csv \
        --method  ${method} \
        --groupby ${groupby} \
        --plots   de_plots/
    """
}
