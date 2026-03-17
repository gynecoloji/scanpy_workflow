"""
Gene set loading utilities.

Supports four sources:
  - "progeny"          : 14 cancer signalling pathways via decoupler.get_progeny()
  - "msigdb_hallmark"  : MSigDB Hallmark (H) collection
  - "msigdb_kegg"      : MSigDB KEGG pathway collection
  - "msigdb_reactome"  : MSigDB Reactome pathway collection
  - "msigdb_gobp"      : MSigDB GO Biological Process collection
  - "custom"           : user-supplied GMT or CSV/TSV file

All loaders return a DataFrame with at minimum:
  source  (str) — gene set name
  target  (str) — gene symbol
  weight  (float, optional) — used by mlm / ulm
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


# Mapping from short source key to MSigDB collection name used by decoupler
_MSIGDB_COLLECTIONS: dict[str, str] = {
    "msigdb_hallmark":  "hallmark",
    "msigdb_kegg":      "kegg_pathways",
    "msigdb_reactome":  "reactome_pathways",
    "msigdb_gobp":      "go_biological_process",
}


def load_genesets(
    source: str,
    organism: str = "human",
    custom_path: str | None = None,
    progeny_top: int = 100,
) -> pd.DataFrame:
    """Return a tidy gene-set DataFrame (source, target[, weight]).

    Parameters
    ----------
    source:
        One of "progeny", "msigdb_hallmark", "msigdb_kegg", "msigdb_reactome",
        "msigdb_gobp", or "custom".
    organism:
        "human" or "mouse" — forwarded to decoupler resource fetchers.
    custom_path:
        Path to a .gmt or .csv/.tsv file.  Required when source == "custom".
    progeny_top:
        Number of top genes per pathway to use for PROGENy.
    """
    if source == "custom":
        if custom_path is None:
            raise ValueError("custom_path is required when source='custom'")
        return _load_custom_file(custom_path)

    if source == "progeny":
        return _load_progeny(organism, progeny_top)

    if source in _MSIGDB_COLLECTIONS:
        return _load_msigdb(source, organism)

    raise ValueError(
        f"Unknown gene-set source: '{source}'. "
        f"Valid values: progeny, {', '.join(_MSIGDB_COLLECTIONS)}, custom"
    )


# ── private loaders ────────────────────────────────────────────────────────────

def _load_progeny(organism: str, top: int) -> pd.DataFrame:
    try:
        import decoupler as dc
    except ImportError:
        raise ImportError(
            "decoupler-py is required for source='progeny'. "
            "Install with: pip install decoupler-py"
        )
    net = dc.get_progeny(organism=organism, top=top)
    # decoupler returns: source, target, weight
    return net[["source", "target", "weight"]].copy()


def _load_msigdb(source: str, organism: str) -> pd.DataFrame:
    try:
        import decoupler as dc
    except ImportError:
        raise ImportError(
            "decoupler-py is required for MSigDB sources. "
            "Install with: pip install decoupler-py"
        )
    collection = _MSIGDB_COLLECTIONS[source]
    msigdb = dc.get_resource("MSigDB", organism=organism)
    subset = msigdb[msigdb["collection"] == collection].copy()
    if subset.empty:
        raise RuntimeError(
            f"MSigDB collection '{collection}' returned no gene sets for organism='{organism}'."
        )
    subset = subset.rename(columns={"geneset": "source", "genesymbol": "target"})
    return subset[["source", "target"]].copy()


def _load_custom_file(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Gene-set file not found: {path}")

    if p.suffix.lower() == ".gmt":
        return _parse_gmt(p)

    # CSV / TSV
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep)

    # Accept either (gene_set / geneset / source) and (gene / genesymbol / target)
    col_map: dict[str, str] = {}
    for cand in ["gene_set", "geneset", "source"]:
        if cand in df.columns:
            col_map[cand] = "source"
            break
    for cand in ["gene", "genesymbol", "target"]:
        if cand in df.columns:
            col_map[cand] = "target"
            break

    if "source" not in col_map.values() or "target" not in col_map.values():
        raise ValueError(
            f"Cannot identify gene-set and gene columns in {path}. "
            "Expected columns: gene_set/geneset/source AND gene/genesymbol/target"
        )
    df = df.rename(columns=col_map)

    cols = ["source", "target"]
    if "weight" in df.columns:
        cols.append("weight")
    return df[cols].copy()


def _parse_gmt(path: Path) -> pd.DataFrame:
    """Parse MSigDB-style GMT (tab-separated: name, description, gene, gene, ...)."""
    rows: list[dict] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gs_name = parts[0]
            genes = [g for g in parts[2:] if g]  # skip description at index 1
            for gene in genes:
                rows.append({"source": gs_name, "target": gene})
    if not rows:
        raise ValueError(f"No gene sets parsed from {path}")
    return pd.DataFrame(rows)
