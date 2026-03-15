// workflow/modules/pca.nf
process PCA {
    label     'scanpy'
    publishDir "${params.output_dir}/pca", mode: 'copy'

    input:
    path h5ad
    val  n_comps
    val  scvi_path   // "true" | "false" — skips scale, asserts adata.X == log_norm

    output:
    path "pca.h5ad"

    script:
    """
    scanpy-workflow pca \
        --input     ${h5ad} \
        --output    pca.h5ad \
        --n-comps   ${n_comps} \
        --scvi-path ${scvi_path}
    """
}
