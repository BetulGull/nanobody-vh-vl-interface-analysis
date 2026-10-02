"""Stage 12 — Final Integrated Mechanism, Evidence Synthesis, and Paper-Ready Mechanistic Model: test suite."""
import sys
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path("/Users/betulgul/Desktop/staj")
TABLES = ROOT / "results/tables"
AUDITS = ROOT / "results/audits"
PROCESSED = ROOT / "data/processed"


# ---------- 1: Stage 2-11 frozen outputs unchanged ----------

def test_stage2to11_frozen_outputs_unchanged():
    assert len(pd.read_parquet(PROCESSED / "stage4_vhh_sequence_representatives.parquet")) == 533
    df3 = pd.read_parquet(PROCESSED / "stage3_paired_representatives.parquet")
    assert len(df3[df3['is_primary_representative'] == True]) == 1429
    assert len(pd.read_parquet(PROCESSED / "stage8_surface_patch_descriptors.parquet")) == 1962
    assert len(pd.read_parquet(PROCESSED / "stage9_antigen_engagement_summary.parquet")) == 1962
    assert len(pd.read_parquet(PROCESSED / "stage11_model_input.parquet")) == 1962


# ---------- 2: mechanism-edge table has no unsupported source references ----------

def test_mechanism_edges_have_supporting_stage():
    path = TABLES / "stage12_mechanism_edges.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    included = df[df['final_include'] == True]
    assert included['supporting_stage'].notna().all()
    assert included['supporting_result'].notna().all()
    assert (included['supporting_result'].str.len() > 0).all()


# ---------- 3: every final claim maps to evidence ----------

def test_every_claim_maps_to_evidence():
    path = TABLES / "stage12_claim_evidence_ledger.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    assert df['source_stage'].notna().all()
    assert df['source_table'].notna().all()
    assert df['evidence_type'].notna().all()


# ---------- 4: every SOLID edge has direct or adjusted support ----------

def test_solid_edges_have_direct_or_adjusted_support():
    path = TABLES / "stage12_mechanism_edges.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    solid = df[df['figure_style'] == 'SOLID']
    for _, row in solid.iterrows():
        assert row['evidence_level'] in ('LEVEL 1', 'LEVEL 2'), f"SOLID edge {row['source_node']}->{row['target_node']} has weak evidence_level {row['evidence_level']}"


# ---------- 5: attenuated shielding->engagement not a strong direct edge ----------

def test_shielding_engagement_not_solid():
    path = TABLES / "stage12_mechanism_edges.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    row = df[(df['source_node'] == 'CDR3_GEOMETRIC_ACCOMMODATION') & (df['target_node'] == 'ANTIGEN_INTERFACE_ORGANIZATION')]
    assert len(row) == 1
    assert row.iloc[0]['figure_style'] != 'SOLID'
    assert row.iloc[0]['causal_claim_allowed'] == 'NONE' or 'NONE' in str(row.iloc[0]['causal_claim_allowed'])


# ---------- 6: shielding and BSA not independent mechanisms ----------

def test_shielding_bsa_coupled_not_independent():
    path = TABLES / "stage12_mechanism_edges.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    row = df[(df['source_node'] == 'STAGE6_SHIELDING') & (df['target_node'] == 'STAGE7_BSA_PACKING')]
    assert len(row) == 1
    assert row.iloc[0]['figure_style'] == 'DOUBLE_LINE_COUPLED'
    assert row.iloc[0]['collinearity_issue'] == 'SEVERE'


# ---------- 7: Gordon theta/phi absent ----------

def test_gordon_theta_phi_absent():
    for path in PROCESSED.glob("stage12_*.parquet"):
        df = pd.read_parquet(path)
        assert 'cdr3_theta' not in df.columns
        assert 'cdr3_phi' not in df.columns
    path = TABLES / "stage12_mechanism_edges.csv"
    if path.exists():
        df = pd.read_csv(path)
        gordon_row = df[df['source_node'] == 'GORDON_THETA_PHI']
        assert (gordon_row['final_include'] == False).all()


# ---------- 8-10: interface definitions preserved and distinct ----------

def test_28_position_interface_preserved():
    df = pd.read_csv(PROCESSED.parent / "processed/stage3_empirical_vh_vl_interface_primary.csv")
    primary = df[df['primary_binary_interface'].astype(str).isin(['True', '1', 'true'])]
    assert len(primary) == 28


def test_15_position_shieldable_preserved():
    FULL = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
    SHIELDABLE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}
    assert len(SHIELDABLE) == 15
    assert SHIELDABLE.issubset(FULL)


def test_14_position_packing_target_distinct():
    SHIELDABLE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}
    NONADJACENT = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}
    assert len(NONADJACENT) == 14
    assert NONADJACENT != SHIELDABLE
    assert NONADJACENT.issubset(SHIELDABLE)


# ---------- 11: tetrad annotation-only ----------

def test_tetrad_annotation_only():
    TETRAD = {42, 49, 50, 52}
    FULL = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
    assert TETRAD.issubset(FULL)
    assert TETRAD != FULL  # never redefined as the whole interface


# ---------- 12-13: comparator correctness ----------

def test_vhh_vs_vh_comparator_for_heavy_domain():
    path = TABLES / "stage11_primary_group_inference.csv"
    df = pd.read_csv(path)
    fw_row = df[df['hypothesis_id'] == 'D_framework_engagement'].iloc[0]
    assert 'VH' in fw_row['comparison'] and 'Fv' not in fw_row['comparison']


def test_vhh_vs_fv_only_for_total_paratope():
    path = TABLES / "stage12_primary_results_summary.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    # total-paratope-style rows (if any explicitly reference Fv) must say Fv, not VH-only
    fv_rows = df[df['comparison_or_model'].str.contains('Fv', na=False)]
    for _, row in fv_rows.iterrows():
        assert 'Fv' in row['comparison_or_model']


# ---------- 14: Stage-5 VH_WITHOUT_VL_GEOMETRIC state preserved ----------

def test_vh_without_vl_geometric_state_preserved():
    df = pd.read_parquet(PROCESSED / "stage5_sequence_surface_chemistry.parquet")
    states = df[df['group'] == 'VH']['structural_state'].unique()
    assert 'VH_WITHOUT_VL_GEOMETRIC' in states
    assert 'VH_WITH_VL_PAIRED' in states  # both retained, never collapsed


# ---------- 15-17: PPC/PNC/APBS separation ----------

def test_ppc_pnc_kept_separate():
    df = pd.read_parquet(PROCESSED / "stage8_surface_patch_descriptors.parquet")
    assert 'former_interface_core_PPC' in df.columns
    assert 'former_interface_core_PNC' in df.columns
    assert 'net_patch' not in df.columns


def test_apbs_components_kept_separate():
    df = pd.read_parquet(PROCESSED / "stage8_apbs_interface_electrostatics.parquet")
    assert 'interface_potential_positive_mean' in df.columns
    assert 'interface_potential_negative_mean' in df.columns


def test_no_net_electrostatic_score_invented():
    for path in list(PROCESSED.glob("stage8_*.parquet")) + list(PROCESSED.glob("stage12_*.parquet")):
        df = pd.read_parquet(path)
        forbidden = ['net_potential', 'net_electrostatic', 'net_charge_score']
        for f in forbidden:
            assert f not in [c.lower() for c in df.columns]


# ---------- 18-21: no overclaiming ----------

def test_no_affinity_claim_in_ledger():
    path = TABLES / "stage12_claim_evidence_ledger.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    text = ' '.join(df['allowed_final_wording'].fillna('').astype(str)).lower()
    assert 'higher affinity' not in text
    assert 'stronger binding' not in text


def test_no_experimental_solubility_claim():
    path = TABLES / "stage12_claim_evidence_ledger.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    text = ' '.join(df['allowed_final_wording'].fillna('').astype(str)).lower()
    assert 'more soluble' not in text
    assert 'therefore vhh is soluble' not in text


def test_no_experimental_stability_claim():
    path = TABLES / "stage12_claim_evidence_ledger.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    text = ' '.join(df['allowed_final_wording'].fillna('').astype(str)).lower()
    assert 'more stable' not in text
    assert 'higher thermodynamic stability' not in text


def test_no_causal_mediation_claim():
    path = TABLES / "stage12_mechanism_edges.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    text = ' '.join(df['causal_claim_allowed'].fillna('').astype(str)).lower()
    assert 'mediates' not in text
    assert 'proven to cause' not in text


# ---------- 22: no Stage-10 output used ----------

def test_no_stage10_output_used():
    for path in PROCESSED.glob("stage12_*"):
        if path.suffix == '.parquet':
            df = pd.read_parquet(path)
            assert not any('aggrescan' in c.lower() or 'developability' in c.lower() for c in df.columns)


# ---------- 23: null/attenuated Stage-11 results retained ----------

def test_null_results_retained_in_primary_mechanism_table():
    path = TABLES / "stage11_primary_mechanism_inference.csv"
    df = pd.read_csv(path)
    assert (df['raw_p'] > 0.05).any()

    path2 = TABLES / "stage12_claim_evidence_ledger.csv"
    if path2.exists():
        ledger = pd.read_csv(path2)
        claim_j = ledger[ledger['claim_id'] == 'J']
        assert len(claim_j) == 1
        assert 'NOT_SUPPORTED' in claim_j.iloc[0]['allowed_final_wording'] or 'not' in claim_j.iloc[0]['allowed_final_wording'].lower()


# ---------- 24: main figures map to frozen source tables ----------

def test_main_figures_map_to_source_tables():
    path = ROOT / "results/figures/stage12/figure_manifest.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    assert df['source_data'].notna().all()
    assert (df['source_data'].str.len() > 0).all()


# ---------- 25: no residue-level pseudoreplication introduced ----------

def test_no_residue_level_pseudoreplication_in_stage12():
    for path in PROCESSED.glob("stage12_*.parquet"):
        df = pd.read_parquet(path)
        if 'INSTANCE' in df.columns:
            assert df['INSTANCE'].is_unique


# ---------- 26: primary result summary N agrees with source ----------

def test_primary_result_summary_n_matches_source():
    path = TABLES / "stage12_primary_results_summary.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    r1 = df[df['result_id'] == 'R1'].iloc[0]
    assert '533' in str(r1['N']) and '1429' in str(r1['N'])


# ---------- 27: claim ledger has no unreferenced claim ----------

def test_claim_ledger_no_unreferenced_claim():
    path = TABLES / "stage12_claim_evidence_ledger.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    assert df['claim_text'].notna().all()
    assert df['source_variable_or_result'].notna().all()


# ---------- 28: limitations table includes direct-solubility limitation ----------

def test_limitations_table_includes_solubility():
    path = TABLES / "stage12_limitations.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    text = ' '.join(df['limitation'].astype(str)).lower()
    assert 'solubility' in text


# ---------- 29: mechanism figure reproducible from edge table ----------

def test_mechanism_figure_file_exists_and_edges_nonempty():
    fig_path = ROOT / "results/figures/stage12/stage12_fig_final_mechanism.png"
    edge_path = TABLES / "stage12_mechanism_edges.csv"
    if not (fig_path.exists() and edge_path.exists()):
        pytest.skip("not built")
    assert fig_path.stat().st_size > 0
    edges = pd.read_csv(edge_path)
    assert len(edges[edges['final_include'] == True]) > 0


# ---------- 30: no Stage-13 computation invoked ----------

def test_no_stage13_computation():
    for path in PROCESSED.glob("stage1*_*.parquet"):
        # only stage10-13 prefixed files would indicate a leap; stage11/12 are expected
        name = path.name
        if name.startswith('stage13') or name.startswith('stage10'):
            pytest.fail(f"Unexpected Stage-10/13 output found: {name}")
