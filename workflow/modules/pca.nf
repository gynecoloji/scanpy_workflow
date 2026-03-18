// workflow/modules/pca.nf
process PCA {
    label     'scanpy'
    publishDir "${params.output_dir}/pca", mode: 'copy'

    input:
    path h5ad
    val  n_comps
    val  scvi_path   // "true" | "false" — skips scale, asserts adata.X == log_norm
    val  evaluate    // "true" | "false"

    output:
    path "pca.h5ad"
    path "evaluation/", optional: true

    script:
    def scvi_flag = (scvi_path == "true") ? "--scvi-path" : ""
    def eval_flag = (evaluate  == "true") ? "--evaluate --eval-dir evaluation/" : ""
    """
    scanpy-workflow pca \
        --input   ${h5ad} \
        --output  pca.h5ad \
        --n-comps ${n_comps} \
        ${scvi_flag} \
        ${eval_flag}
    """
}
