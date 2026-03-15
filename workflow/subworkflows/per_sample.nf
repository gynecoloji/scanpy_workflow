// workflow/subworkflows/per_sample.nf
include { LOAD          } from '../modules/load'
include { AMBIENT       } from '../modules/ambient'
include { QC            } from '../modules/qc'
include { FILTER        } from '../modules/filter'
include { DOUBLETS      } from '../modules/doublets'
include { NORMALIZE       } from '../modules/preprocess'
include { NORMALIZE_SCRAN  } from '../modules/preprocess'
include { HVG_PER_SAMPLE  } from '../modules/preprocess'

workflow PER_SAMPLE {
    take:
    // Channel of [ sample_name, cellranger_dir ]
    samples_ch  // tuple(val, path)

    main:
    // Step 1: Load Cell Ranger output -> .h5ad
    loaded_ch = LOAD(samples_ch)

    // Step 2 (optional): Ambient RNA removal
    if (params.ambient.enabled) {
        // SoupX also needs the raw_feature_bc_matrix dir
        raw_dirs_ch = samples_ch.map { sample, dir ->
            tuple(sample, file("${dir}/raw_feature_bc_matrix"))
        }
        // Join on sample to pair loaded h5ad with its raw dir
        ambient_input_ch = loaded_ch.join(raw_dirs_ch).map { sample, h5ad, raw ->
            tuple(sample, h5ad, raw)
        }
        post_ambient_ch = AMBIENT(
            ambient_input_ch.map { s, h, r -> tuple(s, h) },
            ambient_input_ch.map { s, h, r -> r },
            params.ambient.method
        )
    } else {
        post_ambient_ch = loaded_ch
    }

    // Step 3: QC metrics
    qc_ch = QC(post_ambient_ch)

    // Step 4: Filter
    filtered_ch = FILTER(qc_ch[0])   // qc_ch emits [tuple, path(plots)]

    // Step 5 (optional): Doublet detection
    if (params.doublet.enabled) {
        doublet_ch = DOUBLETS(filtered_ch)
    } else {
        doublet_ch = filtered_ch
    }

    // Step 6: Normalization — branch on method to select correct env
    if (params.normalization.method == "scran") {
        norm_ch = NORMALIZE_SCRAN(doublet_ch)
    } else {
        norm_ch = NORMALIZE(doublet_ch, params.normalization.target_sum)
    }

    // Step 7: Per-sample HVG
    hvg_ch = HVG_PER_SAMPLE(
        norm_ch,
        params.hvg.n_top_genes,
        params.hvg.flavor
    )

    emit:
    preprocessed = hvg_ch   // tuple(sample_name, preprocessed.h5ad)
}
