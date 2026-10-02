"""Stage-0 environment smoke tests for the vhh3d conda environment.

These verify actual imports and basic functional behavior, not just that
conda reports packages as installed.
"""
import io
import yaml
import pytest


def test_core_imports():
    import numpy
    import pandas
    import scipy
    import sklearn
    import statsmodels
    import matplotlib
    import Bio
    import gemmi
    import pyarrow
    import tqdm
    import joblib
    import networkx
    import freesasa  # noqa: F401


def test_gemmi_can_read_real_sabdab_cif():
    import glob
    import gemmi

    candidates = glob.glob("data/raw_extracted/sdh/**/*.cif", recursive=True)
    assert candidates, "no extracted SD-H .cif files found to test gemmi against"

    structure = gemmi.read_structure(candidates[0])
    assert len(structure) >= 1
    assert len(structure[0]) >= 1
    n_atoms = sum(1 for chain in structure[0] for res in chain for atom in res)
    assert n_atoms > 0


def test_freesasa_basic():
    import freesasa

    assert hasattr(freesasa, "Structure")
    assert hasattr(freesasa, "calc")


def test_pandas_parquet_roundtrip(tmp_path):
    import pandas as pd

    df = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})
    path = tmp_path / "test.parquet"
    df.to_parquet(path)
    df2 = pd.read_parquet(path)
    pd.testing.assert_frame_equal(df, df2)


def test_yaml_config_parses():
    sample = """
    seed: 42
    imgt:
      cdr3:
        start: 105
        end: 117
    """
    parsed = yaml.safe_load(io.StringIO(sample))
    assert parsed["seed"] == 42
    assert parsed["imgt"]["cdr3"]["start"] == 105


def test_config_yaml_file_parses():
    from pathlib import Path

    config_path = Path(__file__).resolve().parent.parent / "configs" / "config.yaml"
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    assert cfg["seed"] == 42
    assert cfg["former_interface"]["primary_delta_sasa_threshold_A2"] == 5.0
    assert cfg["former_interface"]["retain_continuous_interface_map"] is True
