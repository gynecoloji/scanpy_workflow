from pathlib import Path
import anndata as ad


def read_h5ad(path: str) -> ad.AnnData:
    """Read an AnnData object from an .h5ad file."""
    return ad.read_h5ad(path)


def write_h5ad(adata: ad.AnnData, path: str) -> None:
    """Write an AnnData object to an .h5ad file, creating parent dirs as needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
