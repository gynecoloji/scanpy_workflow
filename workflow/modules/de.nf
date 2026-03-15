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

process DE_MAST {
    label     'r_env'
    publishDir "${params.output_dir}/de", mode: 'copy'

    input:
    path h5ad
    val  groupby

    output:
    path "de_results.csv"

    script:
    """
    Rscript ${projectDir}/../src/scanpy_workflow/de/mast.R \
        --input   ${h5ad} \
        --groupby ${groupby} \
        --output  de_results.csv
    """
}

process DE_PSEUDOBULK {
    label     'r_env'
    publishDir "${params.output_dir}/de", mode: 'copy'

    input:
    path h5ad
    val  groupby
    val  sample_key

    output:
    path "de_results.csv"

    script:
    """
    Rscript ${projectDir}/../src/scanpy_workflow/de/pseudobulk.R \
        --input      ${h5ad} \
        --groupby    ${groupby} \
        --sample-key ${sample_key} \
        --min-cells  ${params.de.min_cells ?: 10} \
        --output     de_results.csv
    """
}
