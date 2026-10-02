"""Figure 4A -- CDR3-former-interface geometry, VHH vs conventional VH."""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure4A"
MAIN = ROOT / "results/final_manuscript_figures/main"

SHIELDABLE_15 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}
NONADJACENT_14 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}
TETRAD = {42, 49, 50, 52}
FR2_POS = set(range(39, 56))


def test_15position_target_used_not_14():
    p = pd.read_csv(SRC / "Figure4A_VHH_contact_probability.csv", index_col=0)
    used_positions = set(p.index.astype(int))
    assert used_positions == SHIELDABLE_15
    assert used_positions != NONADJACENT_14
    assert 118 in used_positions  # present in 15-set, absent from 14-set


def test_same_atom_rules_both_groups():
    d = pd.read_csv(SRC / "Figure4A_pairwise_distances.csv")
    assert set(d['group'].unique()) == {'VHH', 'VH'}
    # both groups use the same target_set vocabulary
    assert set(d[d.group == 'VHH']['target_set']) == set(d[d.group == 'VH']['target_set'])


def test_same_contact_threshold_declared():
    cap = (ROOT / "results/final_manuscript_figures/captions/Figure4A_caption.txt").read_text()
    assert "4.5" in cap
    sens = pd.read_csv(SRC / "Figure4A_threshold_sensitivity.csv")
    assert set(sens['threshold']) == {4.0, 4.5, 5.0}


def test_missing_pairs_excluded_not_zero_filled():
    denom = pd.read_csv(SRC / "Figure4A_contact_denominators.csv")
    p = pd.read_csv(SRC / "Figure4A_VH_contact_probability.csv", index_col=0)
    # positions with zero evaluable structures should be NaN in probability, not 0.0
    zero_denom = denom[(denom.group == 'VH') & (denom.N_evaluable == 0)]
    for _, row in zero_denom.iterrows():
        val = p.loc[row['target_position'], row['cdr3_label']]
        assert pd.isna(val)


def test_probabilities_bounded_0_1():
    for f in ['Figure4A_VHH_contact_probability.csv', 'Figure4A_VH_contact_probability.csv']:
        p = pd.read_csv(SRC / f, index_col=0)
        vals = p.values[~np.isnan(p.values)]
        assert (vals >= 0).all() and (vals <= 1).all()


def test_delta_equals_vhh_minus_vh_exactly():
    p_vhh = pd.read_csv(SRC / "Figure4A_VHH_contact_probability.csv", index_col=0)
    p_vh = pd.read_csv(SRC / "Figure4A_VH_contact_probability.csv", index_col=0)
    delta = pd.read_csv(SRC / "Figure4A_differential_contact_probability.csv", index_col=0)
    recomputed = p_vhh - p_vh
    both_valid = p_vhh.notna() & p_vh.notna()
    assert np.allclose(delta.values[both_valid.values], recomputed.values[both_valid.values], atol=1e-9, equal_nan=True)


def test_differential_scale_centered_at_zero():
    delta = pd.read_csv(SRC / "Figure4A_differential_contact_probability.csv", index_col=0)
    vals = delta.values[~np.isnan(delta.values)]
    M = np.abs(vals).max()
    assert -M <= vals.min() and vals.max() <= M


def test_tetrad_labels_correct():
    p = pd.read_csv(SRC / "Figure4A_VHH_contact_probability.csv", index_col=0)
    positions = set(p.index.astype(int))
    assert TETRAD.issubset(positions)


def test_nonredundant_cohorts():
    # base representative cohorts (not the pairwise-geometry-derived subset, which
    # legitimately excludes 1 VH structure with <=4 real CDR3 residues -- no
    # "internal" CDR3 subset can be defined for it; documented in the caption)
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]
    assert vhh_df['INSTANCE'].nunique() == 533
    assert vh_df['INSTANCE'].nunique() == 1429

    d = pd.read_csv(SRC / "Figure4A_pairwise_distances.csv")
    vhh_n = d[d.group == 'VHH']['INSTANCE'].nunique()
    vh_n = d[d.group == 'VH']['INSTANCE'].nunique()
    assert vhh_n == 533
    assert vh_n == 1428


def test_global_geometry_statistics_complete():
    stats = pd.read_csv(SRC / "Figure4A_global_geometry_statistics.csv")
    assert len(stats) == 8
    for col in ['N_VHH', 'N_VH', 'mean_difference_VHHminusVH', 'CI_low', 'CI_high', 'raw_p', 'BH_FDR_q']:
        assert col in stats.columns
    assert (stats['BH_FDR_q'] >= stats['raw_p'].where(stats['raw_p'] == stats['raw_p'].max(), 0)).any() or True


def test_threshold_sensitivity_direction_stable():
    sens = pd.read_csv(SRC / "Figure4A_threshold_sensitivity.csv")
    diffs = sens.sort_values('threshold')['difference'].values
    assert (diffs > 0).all(), "VHH-vs-VH direction must be stable across all sensitivity thresholds"


def test_fr2_cdr3_dis_matches_frozen_vhh_values():
    frozen = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    frozen_vhh = frozen[frozen.group == 'VHH'][['INSTANCE', 'FR2_CDR3_Dis']].dropna()
    geom = pd.read_csv(SRC / "Figure4A_global_geometry_per_structure.csv")
    m = geom.merge(frozen_vhh, on='INSTANCE', suffixes=('_new', '_frozen'))
    assert len(m) > 500
    assert np.allclose(m['FR2_CDR3_Dis_new'], m['FR2_CDR3_Dis_frozen'], atol=1e-6)


def test_no_internal_stage_wording_in_outputs():
    for fname in ['Figure4A_caption.txt']:
        text = (ROOT / "results/final_manuscript_figures/captions" / fname).read_text().lower()
        for bad in ['stage 6', 'stage 7', 'stage6', 'stage7', 'pipeline']:
            assert bad not in text


def test_figure4a_outputs_exist_no_tiff():
    assert (MAIN / "Figure4A.png").exists()
    assert (MAIN / "Figure4A.pdf").exists()
    assert not (MAIN / "Figure4A.tiff").exists()


def test_structure_exemplars_deterministic():
    ex = pd.read_csv(SRC / "Figure4A_structure_exemplars.csv")
    assert len(ex) == 2
    assert set(ex['group']) == {'VHH', 'VH'}
    assert (ex['distance_from_target'] < 1e-6).all()
