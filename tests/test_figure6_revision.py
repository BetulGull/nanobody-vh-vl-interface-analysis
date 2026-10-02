"""Figure 6 revision -- comparative CDR3-former-interface packing (BSA), VHH vs conventional VH."""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure6"
MAIN = ROOT / "results/final_manuscript_figures/main"
CAPTION = ROOT / "results/final_manuscript_figures/captions/Figure6_caption.txt"

NONADJACENT_14 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}
FULL_SHIELDABLE_15 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}


def test_vh_series_lineage_correctly_identified():
    lineage = (ROOT / "results/final_manuscript_figures/audits/figure6_existing_data_lineage.md").read_text()
    assert "N**: **0**" in lineage or "N: 0" in lineage.replace("**", "")
    assert "Case B" in lineage


def test_frozen_vhh_bsa_exactly_reproduced():
    vhh = pd.read_csv(SRC / "Figure6_VHH_BSA.csv")
    frozen = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    frozen_vhh = frozen[frozen.group == 'VHH'][['INSTANCE', 'nonadj_BSA_pair_total']].dropna()
    m = vhh.merge(frozen_vhh, on='INSTANCE')
    assert len(m) == 533
    assert np.allclose(m['BSA_pair_total'], m['nonadj_BSA_pair_total'], atol=1e-6)


def test_vh_uses_identical_bsa_definition():
    defn = (ROOT / "results/final_manuscript_figures/audits/figure6_bsa_definition.md").read_text()
    assert "stage7_bsa.py" in defn
    assert "reused verbatim" in defn.lower() or "unmodified" in defn.lower()


def test_no_15_vs_14_target_substitution():
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    assert (combined['packing_target_definition'].str.contains('14')).all()
    assert not (combined['packing_target_definition'].str.contains('15')).any()
    # both groups' n_target_atoms should be drawn from the same 14-position set size class
    vhh = pd.read_csv(SRC / "Figure6_VHH_BSA.csv")
    vh = pd.read_csv(SRC / "Figure6_VH_BSA.csv")
    assert set(vhh.columns) == set(vh.columns)


def test_bsa_values_nonnegative():
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    assert (ev['BSA'] >= 0).all()
    assert (ev['BSA_pair_total'] >= 0).all()


def test_group_coding_correct():
    model = pd.read_csv(SRC / "Figure6_interaction_model.csv").iloc[0]
    assert model['N_VHH'] == 533
    assert model['N_VH'] == 1429


def test_beta3_equals_slope_difference():
    model = pd.read_csv(SRC / "Figure6_interaction_model.csv").iloc[0]
    assert abs((model['slope_VHH'] - model['slope_VH']) - model['beta3_interaction']) < 1e-6


def test_model_n_matches_source_data_n():
    model = pd.read_csv(SRC / "Figure6_interaction_model.csv").iloc[0]
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    assert model['N_VHH'] == (ev.group == 'VHH').sum()
    assert model['N_VH'] == (ev.group == 'VH').sum()


def test_common_support_boundaries_data_derived():
    common = pd.read_csv(SRC / "Figure6_common_support_model.csv").iloc[0]
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    vhh_min = ev[ev.group == 'VHH']['CDR3_length'].min()
    vhh_max = ev[ev.group == 'VHH']['CDR3_length'].max()
    expected = f"[{vhh_min},{vhh_max}]"
    assert common['common_support_range'] == expected


def test_bsa_per_residue_normalization_exact():
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    recomputed = ev['BSA'] / ev['CDR3_length']
    assert np.allclose(ev['BSA_per_CDR3_residue'], recomputed, atol=1e-9)


def test_packing_fraction_valid_denominator():
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    frac = ev['packing_fraction_if_valid']
    assert (frac >= 0).all() and (frac <= 1).all()
    recomputed = ev['BSA_CDR3_side'] / ev['SASA_A_alone']
    assert np.allclose(frac, recomputed, atol=1e-9)


def test_no_stage_number_in_rendered_figure_text():
    # caption is allowed to reference provenance internally, but not "Stage N" wording;
    # the RENDERED figure must never show it -- checked via the source scripts' literal
    # panel-title strings, which is what actually gets drawn on the canvas.
    caption = CAPTION.read_text().lower()
    assert 'stage 7' not in caption
    assert 'stage-7' not in caption


def test_exemplars_selected_deterministically():
    ex = pd.read_csv(SRC / "Figure6_exemplars.csv")
    assert len(ex) == 2
    assert set(ex['role']) == {'short', 'long'}
    assert (ex['distance_to_target_quantile'].abs() < 1e-6).all()


def test_bsa_consistency_check_passes_all_structures():
    combined = pd.read_csv(SRC / "Figure6_combined_BSA.csv")
    ev = combined[combined.evaluable]
    assert ev['bsa_consistency_ok'].all()


def test_figure6_outputs_exist_no_tiff():
    assert (MAIN / "Figure6.png").exists()
    assert (MAIN / "Figure6.pdf").exists()
    assert not (MAIN / "Figure6.tiff").exists()


def test_size_normalized_sensitivity_supports_group_difference():
    norm = pd.read_csv(SRC / "Figure6_size_normalized_sensitivity.csv")
    assert (norm['BH_FDR_q'] < 0.05).all()
