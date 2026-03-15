// workflow/modules/preprocess.nf

// library_size normalization: Python CLI, env_scanpy
process NORMALIZE {
    label     'scanpy'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/preprocessed", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)
    val   target_sum

    output:
    tuple val(sample), path("${sample}_normalized.h5ad")

    script:
    """
    scanpy-workflow normalize \
        --input      ${h5ad} \
        --output     ${sample}_normalized.h5ad \
        --target-sum ${target_sum}
    """
}

// scran normalization: R script, env_r
process NORMALIZE_SCRAN {
    label     'r_env'
    tag       "${sample}"
    publishDir "${params.output_dir}/per_sample/${sample}/preprocessed", mode: 'copy'

    input:
    tuple val(sample), path(h5ad)

    output:
    tuple val(sample), path("${sample}_normalized.h5ad")

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    """
    Rscript ${script_dir}/preprocessing/scran.R \
        --input  ${h5ad} \
        --output ${sample}_normalized.h5ad
    """
}

process HVG_PER_SAMPLE {
    label 'scanpy'
    tag   "${sample}"

    input:
    tuple val(sample), path(h5ad)
    val   n_top_genes
    val   flavor

    output:
    tuple val(sample), path("${sample}_hvg.h5ad")

    script:
    """
    scanpy-workflow hvg \
        --input       ${h5ad} \
        --output      ${sample}_hvg.h5ad \
        --mode        per-sample \
        --n-top-genes ${n_top_genes} \
        --flavor      ${flavor}
    """
}
