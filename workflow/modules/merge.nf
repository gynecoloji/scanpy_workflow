// workflow/modules/merge.nf
process MERGE {
    label     'scanpy'
    publishDir "${params.output_dir}/merged", mode: 'copy'

    input:
    path h5ad_files     // collected list of per-sample .h5ad
    val  sample_names   // space-separated string of sample names
    val  batch_key      // obs column for batch identity

    output:
    path "merged.h5ad"

    script:
    """
    scanpy-workflow merge \
        --inputs    ${h5ad_files} \
        --samples   ${sample_names} \
        --batch-key ${batch_key} \
        --output    merged.h5ad
    """
}
