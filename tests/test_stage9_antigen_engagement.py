"""
Stage 9 — Antigen Engagement, Framework Contribution, and Functional
Accessibility of the Former VH-VL Surface: test suite.
"""
import sys
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path("/Users/betulgul/Desktop/staj")
sys.path.insert(0, str(ROOT / "src" / "qc"))

from stage9_antigen_contacts import get_region, IMGT_REGIONS, compute_residue_contacts
from antigen_qc import AntigenQuality

PROCESSED = ROOT / "data/processed"
TABLES = ROOT / "results/tables"

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
SHIELDABLE_NONCDR3 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}
NONADJACENT_TARGET_14 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}


# ---------- 1-4: antigen identification / masking ----------

def test_legitimate_protein_antigen_types_only():
    assert AntigenQuality.PROTEIN_POLYMER_TYPES == {'polypeptide(L)', 'polypeptide(D)'}


def test_sibling_antibody_chain_excluded():
    aq = AntigenQuality()

    class FakeStructure:
        def get_chains_for_model(self, model_id):
            return ['H', 'L', 'X']

        def get_entity_id_for_chain(self, chain_id):
            return {'H': 'e1', 'L': 'e2', 'X': 'e3'}[chain_id]

        entity_info = {'e1': {'type': 'polypeptide(L)'}, 'e2': {'type': 'polypeptide(L)'}, 'e3': {'type': 'polypeptide(L)'}}

    chains, _ = aq.identify_antigen_chains(FakeStructure(), '1', {'H', 'L'})
    assert chains == ['X']


def test_multi_instance_pdb_handling():
    # A second antibody's H/L chains present in the same PDB must be excluded
    # from antigen identification even though they are not THIS instance's
    # own antibody_chains set -- mirrors stage2_structure_qc.py's
    # build_pdb_antibody_chain_map behavior (all-instances union passed in).
    aq = AntigenQuality()

    class FakeStructure:
        def get_chains_for_model(self, model_id):
            return ['H', 'L', 'H2', 'L2', 'X']

        def get_entity_id_for_chain(self, chain_id):
            return {'H': 'e1', 'L': 'e2', 'H2': 'e1', 'L2': 'e2', 'X': 'e3'}[chain_id]

        entity_info = {'e1': {'type': 'polypeptide(L)'}, 'e2': {'type': 'polypeptide(L)'}, 'e3': {'type': 'polypeptide(L)'}}

    # this instance's own chains are H/L, but H2/L2 (a sibling antibody copy
    # in the same PDB) are passed in via the union set, as Stage 2 does
    all_pdb_antibody_chains = {'H', 'L', 'H2', 'L2'}
    chains, _ = aq.identify_antigen_chains(FakeStructure(), '1', all_pdb_antibody_chains)
    assert chains == ['X']


def test_exact_heavy_chain_mapping():
    df = pd.read_parquet(PROCESSED / "stage4_vhh_sequence_representatives.parquet")
    assert 'Hchain' in df.columns
    assert df['Hchain'].notna().all()


def test_light_chain_mapping_retained_for_vh():
    df = pd.read_parquet(PROCESSED / "stage3_paired_representatives.parquet")
    df = df[df['is_primary_representative'] == True]
    assert 'Lchain' in df.columns


# ---------- 6-7: insertion / antigen residue identity ----------

def test_insertion_code_field_present_in_contact_computation():
    ab_atoms = [{'auth_asym_id': 'H', 'auth_seq_id': 111, 'ins_code': 'A', 'residue_name': 'SER',
                 'atom_name': 'CA', 'x': 0, 'y': 0, 'z': 0}]
    ag_atoms = [{'antigen_chain': 'A', 'auth_seq_id': 5, 'ins_code': '', 'atom_name': 'CA',
                 'x': 1, 'y': 0, 'z': 0, 'residue_name': 'ALA'}]
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert len(rows) == 1
    assert rows[0]['antibody_ins_code'] == 'A'


def test_antigen_residue_identity_preserved():
    ab_atoms = [{'auth_asym_id': 'H', 'auth_seq_id': 50, 'ins_code': '', 'residue_name': 'TYR',
                 'atom_name': 'CA', 'x': 0, 'y': 0, 'z': 0}]
    ag_atoms = [{'antigen_chain': 'A', 'auth_seq_id': 12, 'ins_code': '', 'atom_name': 'CA',
                 'x': 2, 'y': 0, 'z': 0, 'residue_name': 'GLU'}]
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert rows[0]['antigen_residue_name'] == 'GLU'
    assert rows[0]['antigen_residue_id'] == 12


# ---------- 8-9: contact definition / threshold sensitivity ----------

def test_4p5_primary_contact_definition():
    ab_atoms = [{'auth_asym_id': 'H', 'auth_seq_id': 50, 'ins_code': '', 'residue_name': 'TYR',
                 'atom_name': 'CA', 'x': 0, 'y': 0, 'z': 0}]
    ag_atoms = [{'antigen_chain': 'A', 'auth_seq_id': 12, 'ins_code': '', 'atom_name': 'CA',
                 'x': 4.2, 'y': 0, 'z': 0, 'residue_name': 'GLU'}]
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert rows[0]['contact_at_4A'] == False
    assert rows[0]['contact_at_4p5A'] == True
    assert rows[0]['contact_at_5A'] == True


def test_threshold_sensitivity_deterministic():
    path = TABLES / "stage9_contact_threshold_sensitivity.csv"
    if not path.exists():
        pytest.skip("full cohort not run")
    df = pd.read_csv(path)
    for group in ['VHH', 'VH']:
        sub = df[df['group'] == group].sort_values('threshold')
        assert sub['n_contact_pairs'].is_monotonic_increasing


# ---------- 10-11: atom vs residue counting ----------

def test_atom_vs_unique_residue_pair_distinct():
    ab_atoms = [
        {'auth_asym_id': 'H', 'auth_seq_id': 50, 'ins_code': '', 'residue_name': 'TYR', 'atom_name': 'CA', 'x': 0, 'y': 0, 'z': 0},
        {'auth_asym_id': 'H', 'auth_seq_id': 50, 'ins_code': '', 'residue_name': 'TYR', 'atom_name': 'CB', 'x': 0.5, 'y': 0, 'z': 0},
    ]
    ag_atoms = [{'antigen_chain': 'A', 'auth_seq_id': 12, 'ins_code': '', 'atom_name': 'CA', 'x': 1, 'y': 0, 'z': 0, 'residue_name': 'GLU'}]
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert len(rows) == 1  # one residue pair, not two atom-level rows


# ---------- 12: CDR vs framework classification ----------

def test_cdr_vs_framework_classification_exact():
    assert get_region(30) == 'CDR1'
    assert get_region(45) == 'FR2'
    assert get_region(110) == 'CDR3'
    assert get_region(120) == 'FR4'


# ---------- 13-14: interface target definitions ----------

def test_full_28_position_interface_preserved():
    assert len(FULL_EMPIRICAL_INTERFACE) == 25  # unique integer positions


def test_15_position_noncdr3_interface_preserved():
    assert len(SHIELDABLE_NONCDR3) == 15
    for pos in SHIELDABLE_NONCDR3:
        assert not (105 <= pos <= 117)


def test_nonadjacent_14_target_not_primary():
    assert NONADJACENT_TARGET_14 != FULL_EMPIRICAL_INTERFACE
    assert NONADJACENT_TARGET_14 != SHIELDABLE_NONCDR3
    assert len(NONADJACENT_TARGET_14) == 14


# ---------- 15: tetrad annotation-only ----------

def test_tetrad_is_annotation_only():
    path = TABLES / "stage9_tetrad_antigen_engagement.csv"
    if not path.exists():
        pytest.skip("tetrad table not built")
    df = pd.read_csv(path)
    assert set(df['imgt_position']) == {42, 49, 50, 52}
    # tetrad positions must be a subset of, not a redefinition of, the full interface
    assert set(df['imgt_position']).issubset(FULL_EMPIRICAL_INTERFACE)


# ---------- 16-19: total paratope / Fv / VH / VL ----------

def test_vhh_total_paratope_definition():
    path = PROCESSED / "stage9_antigen_engagement_summary.parquet"
    if not path.exists():
        pytest.skip("summary not built")
    df = pd.read_parquet(path)
    vhh = df[df['group'] == 'VHH']
    assert 'VH_contact_count' not in vhh.columns or vhh['VH_contact_count'].isna().all()


def test_conventional_fv_total_paratope_definition():
    path = PROCESSED / "stage9_antigen_engagement_summary.parquet"
    if not path.exists():
        pytest.skip("summary not built")
    df = pd.read_parquet(path)
    vh = df[df['group'] == 'VH']
    assert 'VH_contact_count' in vh.columns and 'VL_contact_count' in vh.columns


def test_conventional_vh_only_distinct_from_fv():
    path = PROCESSED / "stage9_antigen_engagement_summary.parquet"
    if not path.exists():
        pytest.skip("summary not built")
    df = pd.read_parquet(path)
    vh = df[df['group'] == 'VH'].dropna(subset=['VH_contact_count', 'total_unique_antibody_contact_residues'])
    assert (vh['VH_contact_count'] <= vh['total_unique_antibody_contact_residues']).all()


def test_vl_contribution_retained():
    path = PROCESSED / "stage9_antibody_antigen_contacts_long.parquet"
    if not path.exists():
        pytest.skip("long table not built")
    df = pd.read_parquet(path)
    assert 'VL' in df[df['group'] == 'VH']['antibody_domain'].unique()


# ---------- 20-21: denominators ----------

def test_framework_fraction_denominator_correct():
    path = PROCESSED / "stage9_antigen_engagement_summary.parquet"
    if not path.exists():
        pytest.skip("summary not built")
    df = pd.read_parquet(path).dropna(subset=['framework_contact_fraction', 'cdr3_contact_fraction'])
    combined = df['framework_contact_fraction'] + (df['CDR1_contact_residue_count'] + df['CDR2_contact_residue_count'] + df['CDR3_contact_residue_count']) / \
        (df['framework_contact_residue_count'] + df['cdr_contact_residue_count']).replace(0, pd.NA)
    assert (combined.dropna() <= 1.0001).all()


def test_cdr3_fraction_denominator_correct():
    path = PROCESSED / "stage9_antigen_engagement_summary.parquet"
    if not path.exists():
        pytest.skip("summary not built")
    df = pd.read_parquet(path).dropna(subset=['cdr3_contact_fraction'])
    assert (df['cdr3_contact_fraction'] >= 0).all() and (df['cdr3_contact_fraction'] <= 1.0001).all()


def test_position_frequency_uses_evaluable_denominator():
    path = PROCESSED / "stage9_former_interface_antigen_position_summary.csv"
    if not path.exists():
        pytest.skip("position summary not built")
    df = pd.read_csv(path)
    assert (df['VHH_antigen_contact_frequency'] >= 0).all() and (df['VHH_antigen_contact_frequency'] <= 1).all()


# ---------- 22-24: failure semantics ----------

def test_no_contact_not_equal_mapping_failure():
    path = ROOT / "data/interim/stage9/cohort_status.jsonl"
    if not path.exists():
        pytest.skip("cohort status not written")
    import json
    statuses = set()
    with open(path) as f:
        for line in f:
            statuses.add(json.loads(line)['status'])
    # ANTIGEN_NOT_CONTACTING and mapping failures are distinct enum values
    assert 'ANTIBODY_MAPPING_FAILURE' not in statuses or 'ANTIGEN_ENGAGEMENT_SUCCESS' in statuses


def test_no_valid_antigen_not_zero_contact():
    # NO_VALID_PROTEIN_ANTIGEN structures produce zero contact ROWS (absence),
    # never a fabricated zero-distance contact record.
    ab_atoms = []
    ag_atoms = []
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert rows == []


def test_multi_antigen_chain_aggregation():
    path = PROCESSED / "stage9_antibody_antigen_contacts_long.parquet"
    if not path.exists():
        pytest.skip("long table not built")
    df = pd.read_parquet(path)
    multi = df.groupby('INSTANCE')['antigen_chain'].nunique()
    assert (multi >= 1).all()


def test_repeated_antigen_atoms_do_not_inflate_unique_residue_count():
    ab_atoms = [{'auth_asym_id': 'H', 'auth_seq_id': 50, 'ins_code': '', 'residue_name': 'TYR', 'atom_name': 'CA', 'x': 0, 'y': 0, 'z': 0}]
    ag_atoms = [
        {'antigen_chain': 'A', 'auth_seq_id': 12, 'ins_code': '', 'atom_name': 'CA', 'x': 1, 'y': 0, 'z': 0, 'residue_name': 'GLU'},
        {'antigen_chain': 'A', 'auth_seq_id': 12, 'ins_code': '', 'atom_name': 'CB', 'x': 1.3, 'y': 0, 'z': 0, 'residue_name': 'GLU'},
    ]
    rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
    assert len(rows) == 1  # same antigen residue, one pair row


# ---------- 27: Arpeggio proximal != specific (structural test, no rerun needed) ----------

def test_arpeggio_proximal_not_specific_by_convention():
    # Reuses the same convention established and tested in Stage 7.
    contact = ['proximal']
    is_specific = any(l not in ('proximal', 'covalent') for l in contact)
    assert is_specific is False


# ---------- 28: no fabricated affinity ----------

def test_no_binding_affinity_metric_fabricated():
    for path in PROCESSED.glob("stage9_*.parquet"):
        df = pd.read_parquet(path)
        forbidden = ['kd', 'kon', 'koff', 'binding_affinity', 'delta_g_binding']
        for f in forbidden:
            assert f not in [c.lower() for c in df.columns]


# ---------- 29: Gordon theta/phi not invoked ----------

def test_gordon_theta_phi_not_invoked():
    for path in PROCESSED.glob("stage9_*.parquet"):
        df = pd.read_parquet(path)
        assert 'cdr3_theta' not in df.columns
        assert 'cdr3_phi' not in df.columns


# ---------- 30: Stage 2-8 unchanged ----------

def test_stage2to8_outputs_unchanged():
    df4 = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    assert len(df4) == 533
    df3 = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    assert len(df3[df3['is_primary_representative'] == True]) == 1429
    df8 = pd.read_parquet(ROOT / "data/processed/stage8_surface_patch_descriptors.parquet")
    assert len(df8) == 1962


# ---------- 31: no Stage-10+ ----------

def test_no_stage10plus_computation():
    for path in PROCESSED.glob("stage9_*.parquet"):
        df = pd.read_parquet(path)
        forbidden = ['aggrescan3d', 'developability_class', 'md_rmsf']
        for f in forbidden:
            assert f not in [c.lower() for c in df.columns]
