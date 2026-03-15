// workflow/modules/imputation.nf

// MAGIC: Python CLI, env_scanpy
process IMPUTE {
    label     'scanpy'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad
    val  evaluate     // "true" | "false"

    output:
    path "imputed.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow impute \
        --input    ${h5ad} \
        --output   imputed.h5ad \
        --method   magic \
        --evaluate ${evaluate} \
        --eval-dir evaluation/
    """
}

// scVI imputation: Python CLI, env_scvi (requires torch + scvi-tools)
process IMPUTE_SCVI {
    label     'scvi'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad
    val  evaluate

    output:
    path "imputed.h5ad"
    path "evaluation/", optional: true

    script:
    """
    scanpy-workflow impute \
        --input    ${h5ad} \
        --output   imputed.h5ad \
        --method   scvi \
        --evaluate ${evaluate} \
        --eval-dir evaluation/
    """
}

// ALRA imputation: R script, env_r
process IMPUTE_ALRA {
    label     'r_env'
    publishDir "${params.output_dir}/imputation", mode: 'copy'

    input:
    path h5ad

    output:
    path "imputed.h5ad"

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/imputation/alra.R \
        --input  ${h5ad} \
        --output imputed.h5ad
    """
}
