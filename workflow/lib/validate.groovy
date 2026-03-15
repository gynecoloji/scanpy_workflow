// workflow/lib/validate.groovy
// Config validation called from main.nf at workflow startup.

def validateConfig(params) {
    // Hard error: scVI cannot be used for both batch correction and imputation
    boolean bc_scvi  = params.batch_correction.enabled &&
                       params.batch_correction.method == "scvi"
    boolean imp_scvi = params.imputation.enabled &&
                       params.imputation.method == "scvi"
    if (bc_scvi && imp_scvi) {
        error """
        Invalid configuration: scVI cannot be used for both batch correction and imputation simultaneously.
        Set batch_correction.method to 'harmony' or 'bbknn', or set imputation.method to 'magic' or 'alra'.
        """
    }

    // Hard error: scType requires a markers file
    if (params.annotation.automated &&
        params.annotation.method == "sctype" &&
        params.annotation.sctype_markers == null) {
        error "Invalid configuration: annotation.sctype_markers must be set when annotation.method is 'sctype'."
    }

    // Warning: de.groupby column won't exist at startup (created by clustering)
    if (params.de.groupby && params.de.groupby != params.clustering.algorithm) {
        log.warn "de.groupby ('${params.de.groupby}') differs from clustering.algorithm ('${params.clustering.algorithm}'). Ensure the column will be created before DE step."
    }
}
