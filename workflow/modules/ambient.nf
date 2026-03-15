// workflow/modules/ambient.nf
process AMBIENT {
    label 'r_env'
    tag   "${sample}"

    input:
    tuple val(sample), path(h5ad)
    path  raw_dir      // CellRanger raw_feature_bc_matrix/ (only used by SoupX)
    val   method       // "soupx" | "decontx"

    output:
    tuple val(sample), path("${sample}_ambient.h5ad")

    script:
    def script_dir = "${projectDir}/../src/scanpy_workflow"
    if (method == "soupx") {
        """
        Rscript ${script_dir}/ambient/soupx.R \
            --filtered-h5ad ${h5ad} \
            --raw-dir        ${raw_dir} \
            --output         ${sample}_ambient.h5ad
        """
    } else if (method == "decontx") {
        """
        Rscript ${script_dir}/ambient/decontx.R \
            --input  ${h5ad} \
            --output ${sample}_ambient.h5ad
        """
    } else {
        error "ambient.nf: unknown method '${method}'. Expected: soupx | decontx"
    }
}
