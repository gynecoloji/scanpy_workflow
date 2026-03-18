// workflow/modules/batch_correction.nf

// Harmony and BBKNN: Python CLI, env_scanpy
process BATCH_CORRECT {
    label     'scanpy'
    publishDir "${params.output_dir}/batch_correction", mode: 'copy'

    input:
    path h5ad
    val  method      // "harmony" | "bbknn"
    val  batch_key
    val  evaluate    // "true" | "false"

    output:
    path "batch_corrected.h5ad"
    path "evaluation/", optional: true

    script:
    def eval_flag = (evaluate == "true") ? "--evaluate --eval-dir evaluation/" : ""
    """
    scanpy-workflow batch-correct \
        --input     ${h5ad} \
        --output    batch_corrected.h5ad \
        --method    ${method} \
        --batch-key ${batch_key} \
        ${eval_flag}
    """
}

// scVI batch correction: Python CLI, env_scvi (requires torch + scvi-tools)
process BATCH_CORRECT_SCVI {
    label     'scvi'
    publishDir "${params.output_dir}/batch_correction", mode: 'copy'

    input:
    path h5ad
    val  batch_key
    val  evaluate    // "true" | "false"

    output:
    path "batch_corrected.h5ad"
    path "evaluation/", optional: true

    script:
    def eval_flag = (evaluate == "true") ? "--evaluate --eval-dir evaluation/" : ""
    """
    scanpy-workflow batch-correct \
        --input     ${h5ad} \
        --output    batch_corrected.h5ad \
        --method    scvi \
        --batch-key ${batch_key} \
        ${eval_flag}
    """
}
