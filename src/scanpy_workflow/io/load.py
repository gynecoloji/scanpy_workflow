from pathlib import Path
import anndata as ad
import scanpy as sc


def load_cellranger(sample_dir: str) -> ad.AnnData:
    """
    Load Cell Ranger output (MEX or HDF5) from a sample directory.

    Tries the following paths in order:
      1. {sample_dir}/filtered_feature_bc_matrix/   (MEX)
      2. {sample_dir}/raw_feature_bc_matrix/         (MEX)
      3. {sample_dir}/filtered_feature_bc_matrix.h5  (HDF5)
      4. {sample_dir}/raw_feature_bc_matrix.h5        (HDF5)
    """
    base = Path(sample_dir)
    mex_paths = [
        base / "filtered_feature_bc_matrix",
        base / "raw_feature_bc_matrix",
    ]
    h5_paths = [
        base / "filtered_feature_bc_matrix.h5",
        base / "raw_feature_bc_matrix.h5",
    ]

    for mex_path in mex_paths:
        if mex_path.is_dir():
            return sc.read_10x_mtx(str(mex_path), var_names="gene_symbols", cache=False)

    for h5_path in h5_paths:
        if h5_path.is_file():
            return sc.read_10x_h5(str(h5_path))

    raise FileNotFoundError(
        f"No Cell Ranger output found in {sample_dir}. "
        "Expected filtered_feature_bc_matrix/ or filtered_feature_bc_matrix.h5"
    )
