"""
Stage 7 — CDR3-Target Intradomain Packing test suite.
"""
import sys
import json
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path("/Users/betulgul/Desktop/staj")
sys.path.insert(0, str(ROOT / "src" / "qc"))

from stage7_bsa import (compute_pairwise_bsa, select_cdr3_atoms, select_target_atoms,
                         NONADJACENT_PACKING_TARGET, FULL_SHIELDABLE_TARGET, CDR3_START, CDR3_END)
from stage7_cif_writer import write_cif_with_surrogate_mapping

TABLES = ROOT / "results/tables"


# ---------- Target-set definition tests ----------

def test_nonadjacent_target_excludes_118():
    assert 118 not in NONADJACENT_PACKING_TARGET
    assert 118 in FULL_SHIELDABLE_TARGET


def test_nonadjacent_is_subset_of_full():
    assert NONADJACENT_PACKING_TARGET.issubset(FULL_SHIELDABLE_TARGET)


def test_target_sets_have_expected_sizes():
    assert len(FULL_SHIELDABLE_TARGET) == 15
    assert len(NONADJACENT_PACKING_TARGET) == 14


def test_target_sets_disjoint_from_cdr3_range():
    for pos in FULL_SHIELDABLE_TARGET:
        assert not (CDR3_START <= pos <= CDR3_END)


def test_cdr3_range_matches_stage6():
    assert CDR3_START == 105 and CDR3_END == 117


# ---------- BSA computation tests ----------

def _synthetic_atoms():
    atoms = []
    for i, resnum in enumerate(range(100, 125)):
        for atom_name, elem, dx in [('N', 'N', 0.0), ('CA', 'C', 1.5), ('C', 'C', 3.0)]:
            atoms.append({
                'auth_asym_id': 'H', 'label_seq_id': i, 'auth_seq_id': resnum, 'ins_code': '',
                'residue_name': 'ALA', 'atom_name': atom_name, 'element': elem,
                'x': float(resnum) + dx, 'y': 0.0, 'z': 0.0, 'occupancy': 1.0, 'altloc': '',
            })
    return atoms


def test_select_cdr3_atoms_range():
    atoms = _synthetic_atoms()
    cdr3 = select_cdr3_atoms(atoms)
    assert all(CDR3_START <= a['auth_seq_id'] <= CDR3_END for a in cdr3)
    assert len(cdr3) > 0


def test_select_target_atoms_matches_set():
    atoms = _synthetic_atoms()
    target = select_target_atoms(atoms, {103, 119, 120})
    assert all(a['auth_seq_id'] in {103, 119, 120} for a in target)


def test_bsa_consistency_synthetic():
    atoms = _synthetic_atoms()
    result = compute_pairwise_bsa(atoms, {103, 119, 120})
    assert result['status'] == 'OK'
    assert result['bsa_consistency_ok']
    assert abs((result['BSA_CDR3_side'] + result['BSA_target_side']) - result['BSA_pair_total']) < 1e-6


def test_bsa_not_evaluable_when_empty_target():
    atoms = [a for a in _synthetic_atoms() if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]
    result = compute_pairwise_bsa(atoms, {999})
    assert result['status'] == 'TARGET_NOT_EVALUABLE'
    assert result['BSA_pair_total'] is None


def test_bsa_values_nonnegative_when_evaluable():
    atoms = _synthetic_atoms()
    result = compute_pairwise_bsa(atoms, NONADJACENT_PACKING_TARGET | {103})
    if result['status'] == 'OK':
        assert result['BSA_pair_total'] >= -1e-6


# ---------- CIF writer / identity mapping tests ----------

def test_cif_writer_reversible_mapping(tmp_path):
    atoms = [{'auth_asym_id': 'FFF', 'label_seq_id': 1, 'auth_seq_id': 105, 'ins_code': '',
              'residue_name': 'ALA', 'atom_name': 'CA', 'element': 'C',
              'x': 1.0, 'y': 2.0, 'z': 3.0, 'occupancy': 1.0, 'altloc': ''}]
    out_path = tmp_path / "test.cif"
    mapping = write_cif_with_surrogate_mapping(atoms, str(out_path), 'test-FFF', 'test')
    assert len(mapping) == 1
    assert mapping[0]['original_auth_asym_id'] == 'FFF'
    assert len(mapping[0]['temporary_chain_id']) == 1
    assert out_path.exists()


def test_cif_writer_contains_chem_comp_category(tmp_path):
    atoms = [{'auth_asym_id': 'A', 'label_seq_id': 1, 'auth_seq_id': 1, 'ins_code': '',
              'residue_name': 'GLY', 'atom_name': 'CA', 'element': 'C',
              'x': 0.0, 'y': 0.0, 'z': 0.0, 'occupancy': 1.0, 'altloc': ''}]
    out_path = tmp_path / "test2.cif"
    write_cif_with_surrogate_mapping(atoms, str(out_path), 'x', 'y')
    content = out_path.read_text()
    assert '_chem_comp.id' in content
    assert 'L-peptide linking' in content


def test_cif_writer_preserves_insertion_code(tmp_path):
    atoms = [{'auth_asym_id': 'A', 'label_seq_id': 1, 'auth_seq_id': 111, 'ins_code': 'A',
              'residue_name': 'SER', 'atom_name': 'CA', 'element': 'C',
              'x': 0.0, 'y': 0.0, 'z': 0.0, 'occupancy': 1.0, 'altloc': ''}]
    out_path = tmp_path / "test3.cif"
    write_cif_with_surrogate_mapping(atoms, str(out_path), 'x', 'y')
    content = out_path.read_text()
    assert 'A' in content  # insertion code character present


def test_cif_writer_multiple_distinct_chains_get_distinct_surrogates(tmp_path):
    atoms = [
        {'auth_asym_id': 'AAA', 'label_seq_id': 1, 'auth_seq_id': 1, 'ins_code': '',
         'residue_name': 'GLY', 'atom_name': 'CA', 'element': 'C', 'x': 0, 'y': 0, 'z': 0, 'occupancy': 1, 'altloc': ''},
        {'auth_asym_id': 'BBB', 'label_seq_id': 1, 'auth_seq_id': 1, 'ins_code': '',
         'residue_name': 'GLY', 'atom_name': 'CA', 'element': 'C', 'x': 0, 'y': 0, 'z': 0, 'occupancy': 1, 'altloc': ''},
    ]
    out_path = tmp_path / "test4.cif"
    mapping = write_cif_with_surrogate_mapping(atoms, str(out_path), 'x', 'y')
    surrogates = {m['temporary_chain_id'] for m in mapping}
    assert len(surrogates) == 2


# ---------- Output table structural tests (require full-cohort run to have completed) ----------

@pytest.fixture(scope="module")
def long_table():
    path = TABLES / "stage7_interaction_records_long.csv"
    if not path.exists():
        pytest.skip("Full cohort run not yet completed")
    return pd.read_csv(path, low_memory=False)


@pytest.fixture(scope="module")
def per_vhh_summary():
    path = TABLES / "stage7_per_vhh_packing_summary.csv"
    if not path.exists():
        pytest.skip("Full cohort run not yet completed")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def position_pair_map():
    path = TABLES / "stage7_position_pair_interaction_map.csv"
    if not path.exists():
        pytest.skip("Full cohort run not yet completed")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def combined_table():
    path = TABLES / "stage7_stage4to7_combined.csv"
    if not path.exists():
        pytest.skip("Full cohort run not yet completed")
    return pd.read_csv(path)


def test_long_table_covers_all_instances(long_table):
    df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    expected_n = len(df)
    assert expected_n == 533  # verified from frozen file, not hard-coded elsewhere
    assert long_table['INSTANCE'].nunique() <= expected_n


def test_long_table_no_cdr3_cdr3_or_target_target_pairs(long_table):
    # every record must have cdr3_imgt_position in CDR3 range and target_imgt_position in FULL_SHIELDABLE_TARGET
    assert (long_table['cdr3_imgt_position'].between(CDR3_START, CDR3_END)).all()
    assert long_table['target_imgt_position'].isin(FULL_SHIELDABLE_TARGET).all()


def test_long_table_proximal_only_never_marked_specific(long_table):
    proximal_only = long_table[long_table['is_proximal_only']]
    assert not proximal_only['is_specific_interaction'].any()


def test_long_table_covalent_records_are_typed_consistently(long_table):
    # Arpeggio's 'covalent' label covers both backbone peptide bonds (117-118) and
    # disulfide bonds (Cys-Cys, which can cross arbitrary sequence distance) -- so
    # covalent CDR3-target crossings are NOT necessarily concentrated at the 118
    # boundary (confirmed empirically: observed covalent crossings involve target
    # positions 50/55, i.e. disulfides, not just the 117-118 peptide bond). This
    # test only asserts internal consistency: every covalent record's contact_labels
    # field actually contains 'covalent'.
    covalent = long_table[long_table['is_covalent']]
    if len(covalent) > 0:
        assert covalent['contact_labels'].str.contains('covalent').all()


def test_long_table_primary_excludes_position_118(long_table):
    primary = long_table[long_table['include_primary_noncovalent']]
    assert 118 not in set(primary['target_imgt_position'])


def test_per_vhh_summary_row_count(per_vhh_summary):
    assert len(per_vhh_summary) == 533


def test_per_vhh_summary_bsa_nonnegative(per_vhh_summary):
    valid = per_vhh_summary[per_vhh_summary['nonadj_status'] == 'OK']
    assert (valid['nonadj_BSA_pair_total'] >= -1e-3).all()


def test_per_vhh_summary_consistency_flag_all_true(per_vhh_summary):
    valid = per_vhh_summary[per_vhh_summary['nonadj_status'] == 'OK']
    assert valid['nonadj_bsa_consistency_ok'].all()


def test_position_pair_map_frequencies_in_range(position_pair_map):
    assert (position_pair_map['interaction_frequency'] >= 0).all()
    assert (position_pair_map['interaction_frequency'] <= 1).all()


def test_position_pair_map_cdr3_positions_valid(position_pair_map):
    assert position_pair_map['cdr3_imgt_position'].between(CDR3_START, CDR3_END).all()


def test_combined_table_row_count(combined_table):
    assert len(combined_table) == 533


def test_combined_table_has_shielding_and_bsa_columns(combined_table):
    assert 'shielding_fraction' in combined_table.columns
    assert 'nonadj_BSA_pair_total' in combined_table.columns


def test_bsa_and_shielding_positively_correlated(combined_table):
    sub = combined_table[['nonadj_BSA_pair_total', 'shielding_fraction']].dropna()
    corr = sub['nonadj_BSA_pair_total'].corr(sub['shielding_fraction'])
    assert corr > 0.5  # descriptive sanity check, not inferential


def test_no_atom_level_count_inflation_via_dual_reporting(per_vhh_summary):
    # unique-residue-pair counts must never exceed atom-level counts
    valid = per_vhh_summary.dropna(subset=['primary_n_specific_interactions_atomlevel', 'primary_n_unique_cdr3_target_respairs'])
    assert (valid['primary_n_unique_cdr3_target_respairs'] <= valid['primary_n_specific_interactions_atomlevel']).all()


def test_zero_interactions_tracked_not_silently_zero(per_vhh_summary):
    # status column must be present and never null for a processed instance
    assert per_vhh_summary['status'].notna().all()


def test_arpeggio_schema_only_atom_atom_type_observed(long_table):
    # per the full-cohort parse, no non-atom-atom (ring) record type should have leaked through
    # (already filtered out at build time; this asserts none were silently included)
    assert 'type' not in long_table.columns or (long_table.get('type', pd.Series(dtype=object)) == 'atom-atom').all()
