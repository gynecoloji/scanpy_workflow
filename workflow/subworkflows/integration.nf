// workflow/subworkflows/integration.nf
include { MERGE              } from '../modules/merge'
include { HVG_POST_MERGE     } from '../modules/hvg_postmerge'
include { SCALE              } from '../modules/scale'
include { PCA                } from '../modules/pca'
include { BATCH_CORRECT      } from '../modules/batch_correction'
include { BATCH_CORRECT_SCVI } from '../modules/batch_correction'
include { IMPUTE             } from '../modules/imputation'
include { IMPUTE_SCVI        } from '../modules/imputation'
include { IMPUTE_ALRA        } from '../modules/imputation'

workflow INTEGRATION {
    take:
    preprocessed_ch   // Channel of tuple(sample_name, preprocessed.h5ad)

    main:
    // Step 8: Merge samples
    // Use multiMap to split channel once, avoiding double-consumption.
    preprocessed_ch
        .multiMap { s, h ->
            h5ads:   h
            samples: s
        }
        .set { split_ch }

    all_h5ads_ch    = split_ch.h5ads.collect()
    all_samples_str = split_ch.samples.collect().map { names -> names.join(' ') }

    merged_h5ad = MERGE(
        all_h5ads_ch,
        all_samples_str,
        params.batch_correction.batch_key
    )

    // Step 9: Post-merge HVG
    hvg_h5ad = HVG_POST_MERGE(
        merged_h5ad,
        params.hvg.post_merge_n_top_genes,
        params.hvg.flavor,
        params.batch_correction.batch_key
    )

    // Step 10: Scale (skip if scVI batch correction — scVI path keeps adata.X = log_norm)
    boolean scvi_bc = (params.batch_correction.enabled &&
                       params.batch_correction.method == "scvi")
    if (scvi_bc) {
        pre_pca_h5ad = hvg_h5ad
    } else {
        pre_pca_h5ad = SCALE(hvg_h5ad)
    }

    // Step 11: PCA
    pca_h5ad = PCA(pre_pca_h5ad, params.pca.n_comps, scvi_bc.toString())

    // Step 12 (optional): Batch correction — branch on method to select env
    if (params.batch_correction.enabled) {
        if (params.batch_correction.method == "scvi") {
            bc_result    = BATCH_CORRECT_SCVI(pca_h5ad, params.batch_correction.batch_key)
        } else {
            bc_result    = BATCH_CORRECT(pca_h5ad, params.batch_correction.method,
                                         params.batch_correction.batch_key)
        }
        post_bc_h5ad = bc_result[0]
    } else {
        post_bc_h5ad = pca_h5ad
    }

    // Step 13 (optional): Imputation — branch on method to select env
    if (params.imputation.enabled) {
        def eval_str = params.imputation.evaluate.toString()
        if (params.imputation.method == "scvi") {
            imp_result = IMPUTE_SCVI(post_bc_h5ad, eval_str)
        } else if (params.imputation.method == "alra") {
            imp_result = IMPUTE_ALRA(post_bc_h5ad)
        } else {
            imp_result = IMPUTE(post_bc_h5ad, eval_str)
        }
        final_h5ad = imp_result[0]
    } else {
        final_h5ad = post_bc_h5ad
    }

    emit:
    integrated = final_h5ad   // path to integrated .h5ad
}
