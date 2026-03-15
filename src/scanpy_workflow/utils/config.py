import warnings
import yaml


def load_config(config_path: str) -> dict:
    """Load params.yaml and validate. Raises ValueError for hard errors."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}
    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    bc = cfg.get("batch_correction", {})
    imp = cfg.get("imputation", {})
    bc_enabled = bc.get("enabled", False)
    imp_enabled = imp.get("enabled", False)
    bc_method = bc.get("method", "harmony")
    imp_method = imp.get("method", "magic")

    if bc_enabled and imp_enabled and bc_method == "scvi" and imp_method == "scvi":
        raise ValueError(
            "Invalid config: batch_correction.method=scvi and imputation.method=scvi "
            "cannot both be enabled. scVI cannot be used for both steps simultaneously."
        )

    de_groupby = cfg.get("de", {}).get("groupby")
    clustering_algo = cfg.get("clustering", {}).get("algorithm", "leiden")
    if de_groupby and de_groupby != clustering_algo:
        warnings.warn(
            f"de.groupby='{de_groupby}' does not match clustering.algorithm='{clustering_algo}'. "
            "Ensure this .obs column will exist at DE runtime.",
            UserWarning,
            stacklevel=3,
        )
