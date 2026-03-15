// workflow/modules/clustering.nf
process CLUSTER {
    label     'scanpy'
    publishDir "${params.output_dir}/clustering", mode: 'copy'

    input:
    path h5ad
    val  algorithm   // "leiden" | "louvain"
    val  resolution
    val  n_neighbors
    val  evaluate    // "true" | "false"

    output:
    path "clustered.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow cluster \
        --input       ${h5ad} \
        --output      clustered.h5ad \
        --algorithm   ${algorithm} \
        --resolution  ${resolution} \
        --n-neighbors ${n_neighbors} \
        --evaluate    ${evaluate} \
        --eval-dir    evaluation/
    """
}
