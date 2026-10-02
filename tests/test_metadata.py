"""Stage 1 metadata-pipeline unit tests."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "qc"))
import stage1_pipeline as sp  # noqa: E402

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")


def test_normalize_method_variants():
    assert sp.normalize_method("X-RAY") == "XRAY"
    assert sp.normalize_method("XRAY") == "XRAY"
    assert sp.normalize_method("x-ray diffraction") == "XRAY"
    assert sp.normalize_method("NA") is None
    assert sp.normalize_method("") is None
    assert sp.normalize_method("EM") == "EM"


def test_coerce_resolution():
    assert sp.coerce_resolution("2.8") == 2.8
    assert sp.coerce_resolution("NA") is None
    assert sp.coerce_resolution("") is None
    assert sp.coerce_resolution("not_a_number") is None


def test_resolution_le_3_logic():
    df = pd.DataFrame({
        "method": ["XRAY", "XRAY", "XRAY", "XRAY"],
        "resolution": ["2.8", "3.0", "3.1", "NA"],
        "antigen_type": ["PROTEIN"] * 4,
        "type": ["SD-H"] * 4,
        "source_cohort": ["SD-H"] * 4,
        "Hchain": ["A", "A", "A", "A"],
        "Lchain": ["NA"] * 4,
    })
    df = sp.add_qc_flags(df)
    assert list(df["qc_resolution_le_3"]) == [True, True, False, False]
    assert df["qc_resolution_numeric"].isna().sum() == 1


def test_camelid_species_classification():
    row = {
        "heavy_taxid": "9838", "heavy_species": "camelus dromedarius",
        "organism_taxid": "9838", "organism": "camelus dromedarius",
        "heavy_subclass": "NA", "compound": "a VHH nanobody complex", "short_header": "IMMUNE SYSTEM",
    }
    result = sp.classify_sdh_row(row)
    assert result["camelid_species_flag"] is True
    assert result["sdh_classification"] == "natural_camelid_candidate"
    assert result["metadata_conflict_flag"] is False


def test_human_sdh_classification_overrides_text():
    row = {
        "heavy_taxid": "9606", "heavy_species": "homo sapiens",
        "organism_taxid": "9606", "organism": "homo sapiens",
        "heavy_subclass": "IGHV3 (Homsap)",
        "compound": "a single domain antibody nanobody-like binder",  # text says nanobody
        "short_header": "IMMUNE SYSTEM",
    }
    result = sp.classify_sdh_row(row)
    # explicit Homo sapiens heavy-chain provenance must win over descriptive nanobody/VHH text
    assert result["sdh_classification"] == "human_sdVH_candidate"
    assert result["camelid_species_flag"] is False


def test_synthetic_sdh_classification():
    row = {
        "heavy_taxid": "32630", "heavy_species": "synthetic construct",
        "organism_taxid": "32630", "organism": "synthetic construct",
        "heavy_subclass": "NA", "compound": "a de novo-designed VHH", "short_header": "DE NOVO PROTEIN",
    }
    result = sp.classify_sdh_row(row)
    assert result["sdh_classification"] == "synthetic_or_de_novo_candidate"


def test_mutant_alone_does_not_exclude_natural_camelid():
    row = {
        "heavy_taxid": "9844", "heavy_species": "lama glama",
        "organism_taxid": "9844", "organism": "lama glama",
        "heavy_subclass": "NA", "compound": "a nanobody point mutant complex", "short_header": "IMMUNE SYSTEM",
    }
    result = sp.classify_sdh_row(row)
    assert result["mutation_flag"] is True
    assert result["sdh_classification"] == "natural_camelid_candidate"
    assert result["metadata_conflict_flag"] is False


def test_construct_normalization_and_scfv_exclusion():
    row_conventional = {"chainsharing_construct": "NA"}
    row_scfv = {"chainsharing_construct": "SCFV"}
    c1 = sp.classify_construct(row_conventional)
    c2 = sp.classify_construct(row_scfv)
    assert c1["construct_conventional"] is True
    assert c2["construct_scfv"] is True
    assert c1["construct_bucket"] == "conventional"
    assert c2["construct_bucket"] == "SCFV"


def test_paired_h_l_metadata_check():
    df = pd.DataFrame({
        "source_cohort": ["FAB", "FAB", "FV"],
        "Hchain": ["H", "H", "H"],
        "Lchain": ["L", "NA", "L"],
        "method": ["XRAY", "XRAY", "XRAY"],
        "resolution": ["2.0", "2.0", "2.0"],
        "antigen_type": ["PROTEIN"] * 3,
        "type": ["FAB", "FAB", "FV"],
        "chainsharing_construct": ["NA", "NA", "NA"],
    })
    df = sp.add_qc_flags(df)
    df = sp.apply_paired_classification(df)
    assert list(df["paired_metadata_candidate"]) == [True, False, True]


def test_protein_antigen_annotation_detection():
    df = pd.DataFrame({
        "method": ["XRAY", "XRAY"],
        "resolution": ["2.0", "2.0"],
        "antigen_type": ["PROTEIN", "PROTEIN|HAPTEN"],
        "type": ["FAB", "FAB"],
        "source_cohort": ["FAB", "FAB"],
        "Hchain": ["H", "H"],
        "Lchain": ["L", "L"],
    })
    df = sp.add_qc_flags(df)
    assert df["qc_contains_protein_annotation"].all()


def test_duplicate_archive_does_not_duplicate_master_rows():
    master_path = PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet"
    if not master_path.exists():
        pytest.skip("master_stage1.parquet not built yet")
    df = pd.read_parquet(master_path)
    n_fab = len(df[df["source_cohort"] == "FAB"])
    assert n_fab == 2798, f"expected 2798 FAB rows, got {n_fab} (duplicate archive must not inflate count)"
    n_total = len(df)
    assert n_total == 4931, f"expected 4931 total rows, got {n_total}"


def test_master_stage1_total_matches_expected():
    master_path = PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet"
    if not master_path.exists():
        pytest.skip("master_stage1.parquet not built yet")
    df = pd.read_parquet(master_path)
    assert len(df) == 1597 + 2798 + 536
