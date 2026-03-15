// workflow/modules/hvg_postmerge.nf
process HVG_POST_MERGE {
    label     'scanpy'
    publishDir "${params.output_dir}/merged", mode: 'copy'

    input:
    path h5ad
    val  n_top_genes
    val  flavor
    val  batch_key     // passed to batch_key param for scib-style HVG

    output:
    path "merged_hvg.h5ad"

    script:
    """
    scanpy-workflow hvg \
        --input        ${h5ad} \
        --output       merged_hvg.h5ad \
        --mode         post-merge \
        --n-top-genes  ${n_top_genes} \
        --flavor       ${flavor} \
        --batch-key    ${batch_key}
    """
}
