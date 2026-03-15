import pytest
import warnings
from scanpy_workflow.utils.config import load_config


def test_load_minimal_config(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: harmony
  batch_key: sample
imputation:
  enabled: false
  method: magic
clustering:
  algorithm: leiden
de:
  method: wilcoxon
  groupby: leiden
""")
    cfg = load_config(str(cfg_file))
    assert cfg["batch_correction"]["method"] == "harmony"


def test_scvi_conflict_raises(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: scvi
imputation:
  enabled: true
  method: scvi
""")
    with pytest.raises(ValueError, match="scVI cannot be used for both"):
        load_config(str(cfg_file))


def test_de_groupby_mismatch_warns(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
clustering:
  algorithm: leiden
de:
  method: wilcoxon
  groupby: louvain
""")
    with pytest.warns(UserWarning, match="does not match clustering.algorithm"):
        load_config(str(cfg_file))


def test_no_scvi_conflict_when_one_disabled(tmp_path):
    cfg_file = tmp_path / "params.yaml"
    cfg_file.write_text("""
batch_correction:
  enabled: true
  method: scvi
imputation:
  enabled: false
  method: scvi
""")
    cfg = load_config(str(cfg_file))  # should not raise
    assert cfg["batch_correction"]["method"] == "scvi"
