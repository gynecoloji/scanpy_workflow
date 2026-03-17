// workflow/modules/pathway.nf
process PATHWAY_SCORE {
    label     'scanpy'
    publishDir "${params.output_dir}/pathway", mode: 'copy'

    input:
    path h5ad
    val  source    // "progeny" | "msigdb_hallmark" | "msigdb_kegg" | "msigdb_reactome" | "msigdb_gobp" | "custom"
    val  method    // "aucell" | "ulm" | "scanpy_score"
    val  groupby
    path custom_genesets   // pass file("NO_FILE") when source != "custom"

    output:
    path "pathway_scored.h5ad"
    path "pathway_results/"

    script:
    def custom_arg = (source == "custom") ? "--custom-genesets ${custom_genesets}" : ""
    """
    scanpy-workflow pathway-score \
        --input      ${h5ad} \
        --output     pathway_scored.h5ad \
        --source     ${source} \
        --method     ${method} \
        --groupby    ${groupby} \
        --min-n      ${params.pathway.min_n ?: 5} \
        --organism   ${params.pathway.organism ?: "human"} \
        --output-dir pathway_results/ \
        ${custom_arg}
    """
}
