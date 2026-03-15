#!/usr/bin/env nextflow
// workflow/main.nf — Scanpy end-to-end single-cell workflow
// DSL2 is the default in Nextflow >= 22.03

include { PER_SAMPLE  } from './subworkflows/per_sample'
include { INTEGRATION } from './subworkflows/integration'
include { CLUSTER     } from './modules/clustering'
include { RANK_GENES  } from './modules/annotation'
include { ANNOTATE_CELLTYPIST } from './modules/annotation'
include { ANNOTATE_SCTYPE     } from './modules/annotation'
include { DE             } from './modules/de'
include { DE_MAST        } from './modules/de'
include { DE_PSEUDOBULK  } from './modules/de'
include { REPORT      } from './modules/report'

// Note: workflow/lib/validate.groovy is auto-loaded by Nextflow; no include needed.

workflow {
    // Config validation (exits with error on invalid combos)
    // validateConfig() is defined in workflow/lib/validate.groovy (auto-loaded)
    validateConfig(params)

    // Build per-sample channel: tuple(sample_name, cellranger_dir)
    samples_ch = Channel.fromList(params.samples)
        .map { sample -> tuple(sample, file("${params.input_dir}/${sample}")) }

    // Steps 1-7: Per-sample processing (parallelized)
    PER_SAMPLE(samples_ch)

    // Steps 8-13: Integration (merge -> HVG -> scale -> PCA -> BC -> imputation)
    INTEGRATION(PER_SAMPLE.out.preprocessed)

    // Steps 14-15: Neighbors + clustering + evaluation
    cluster_result = CLUSTER(
        INTEGRATION.out.integrated,
        params.clustering.algorithm,
        params.clustering.resolution,
        params.clustering.n_neighbors,
        params.clustering.evaluate.toString()
    )
    // Convert to value channel so multiple downstream processes can consume it.
    // Queue channels are single-consumer; .first() re-emits to each subscriber.
    clustered_h5ad = cluster_result[0].first()

    // Step 16a: Unsupervised annotation (rank genes)
    if (params.annotation.unsupervised) {
        rank_result = RANK_GENES(
            clustered_h5ad,
            params.de.groupby ?: params.clustering.algorithm
        )
        annotated_h5ad = rank_result[0]
    } else {
        annotated_h5ad = clustered_h5ad
    }

    // Step 16b: Automated annotation (CellTypist or scType) — runs in parallel with 16a
    // validateConfig() already guarantees sctype_markers is set if method == "sctype"
    if (params.annotation.automated) {
        if (params.annotation.method == "celltypist") {
            ANNOTATE_CELLTYPIST(
                clustered_h5ad,
                params.annotation.celltypist_model,
                params.clustering.algorithm
            )
        } else if (params.annotation.method == "sctype") {
            ANNOTATE_SCTYPE(
                clustered_h5ad,
                file(params.annotation.sctype_markers),
                params.clustering.algorithm
            )
        }
    }

    // Step 17: Differential expression — branch on method to select correct env
    def de_groupby = params.de.groupby ?: params.clustering.algorithm
    if (params.de.method == "mast") {
        de_result = DE_MAST(annotated_h5ad, de_groupby)
    } else if (params.de.method == "pseudobulk") {
        de_result = DE_PSEUDOBULK(
            annotated_h5ad,
            de_groupby,
            params.de.sample_key ?: "sample"
        )
    } else {
        de_result = DE(annotated_h5ad, params.de.method, de_groupby)
    }

    // Step 18: HTML report + final .h5ad
    // de_result[0] is the CSV for all DE methods (mast/pseudobulk emit single output)
    if (params.report) {
        def de_csv = (params.de.method == "mast" || params.de.method == "pseudobulk")
            ? de_result
            : de_result[0]
        REPORT(
            annotated_h5ad,
            de_csv,
            params.output_dir
        )
    }
}
