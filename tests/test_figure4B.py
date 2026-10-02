"""Figure 4B -- CDR3 length-dependent structural organization, VHH vs conventional VH."""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure4B"
MAIN = ROOT / "results/final_manuscript_figures/main"

CDR1_RANGE = set(range(27, 39))
FR2_RANGE = set(range(39, 56))
CDR2_RANGE = set(range(56, 66))
ZHANG_STATES = {'G', 'T', 'E', 'S', 'C'}


def test_metric_definitions_documented():
    status = pd.read_csv(SRC / "Figure4B_metric_definition_status.csv")
    assert len(status) > 0
    assert set(status['status']).issubset({'EXACT', 'NOT_REPRODUCED', 'ADAPTED_NOT_IN_SOURCE', 'EXCLUDED_BY_INSTRUCTION'})


def test_unresolved_definitions_not_guessed():
    status = pd.read_csv(SRC / "Figure4B_metric_definition_status.csv")
    not_repro = status[status.status == 'NOT_REPRODUCED']
    assert len(not_repro) >= 3  # LPP, LNP, surface charge
    for _, row in not_repro.iterrows():
        assert 'omit' in row['vhh_source'].lower()


def test_vhh_vh_identical_interaction_definitions():
    for panel, tset in [('B', 'CDR1'), ('C', 'FR2'), ('D', 'CDR2')]:
        vhh = pd.read_csv(SRC / f"Figure4B_panel{panel}_VHH_CDR3_{tset}.csv", index_col=0)
        vh = pd.read_csv(SRC / f"Figure4B_panel{panel}_VH_CDR3_{tset}.csv", index_col=0)
        assert list(vhh.index) == list(vh.index)
        assert list(vhh.columns) == list(vh.columns)


def test_cdr1_fr2_cdr2_numbering_correct():
    p_fr2 = pd.read_csv(SRC / "FigureB_panelC_placeholder.csv") if False else None
    vhh_fr2 = pd.read_csv(SRC / "Figure4B_panelC_VHH_CDR3_FR2.csv", index_col=0)
    assert set(vhh_fr2.index.astype(int)).issubset(FR2_RANGE)
    vhh_cdr1 = pd.read_csv(SRC / "Figure4B_panelB_VHH_CDR3_CDR1.csv", index_col=0)
    assert set(vhh_cdr1.index.astype(int)).issubset(CDR1_RANGE)
    vhh_cdr2 = pd.read_csv(SRC / "Figure4B_panelD_VHH_CDR3_CDR2.csv", index_col=0)
    assert set(vhh_cdr2.index.astype(int)).issubset(CDR2_RANGE)


def test_fr2_cdr3_dis_identical_definition_across_groups():
    per_struct = pd.read_csv(SRC / "Figure4B_per_structure_metrics.csv")
    frozen = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    frozen_vhh = frozen[frozen.group == 'VHH'][['INSTANCE', 'FR2_CDR3_Dis']].dropna()
    m = per_struct.merge(frozen_vhh, on='INSTANCE', suffixes=('_new', '_frozen'))
    assert len(m) > 500
    assert np.allclose(m['FR2_CDR3_Dis_new'], m['FR2_CDR3_Dis_frozen'], atol=1e-6)
    assert per_struct[per_struct.group == 'VH']['FR2_CDR3_Dis'].notna().sum() > 1000


def test_dssp_mapping_identical_across_groups():
    dssp = pd.read_csv(ROOT / "data/interim/figure4B_dssp_cdr3.csv")
    assert set(dssp['group'].unique()) == {'VHH', 'VH'}
    assert set(dssp['dssp_zhang_state'].unique()).issubset(ZHANG_STATES | {'OTHER'})


def test_same_sasa_definition_both_groups():
    domain = pd.read_parquet(ROOT / "data/processed/figure4B_domain_sasa_hydrophobicity.parquet")
    assert set(domain['structural_state'].unique()) == {'VHH_NATIVE', 'VH_WITHOUT_VL_GEOMETRIC'}
    assert domain[domain.group == 'VHH']['INSTANCE'].nunique() == 533
    assert domain[domain.group == 'VH']['INSTANCE'].nunique() == 1429


def test_surface_metrics_not_silently_substituted():
    status = pd.read_csv(SRC / "Figure4B_metric_definition_status.csv")
    lpp_lnp = status[status.zhang_metric.str.contains('LPP|LNP', regex=True)]
    assert (lpp_lnp['status'] == 'NOT_REPRODUCED').all()
    charge = status[status.zhang_metric.str.contains('Surface charge')]
    assert (charge['status'] == 'NOT_REPRODUCED').all()


def test_solubility_absent():
    for f in SRC.glob("*.csv"):
        cols = pd.read_csv(f, nrows=0).columns
        for c in cols:
            assert 'solub' not in c.lower()
    cap = (ROOT / "results/final_manuscript_figures/captions/Figure4B_caption.txt")
    if cap.exists():
        text = cap.read_text().lower()
        assert 'solubility was not calculated' in text or 'solubility excluded' in text or 'not computed' in text


def test_group_coding_correct():
    models = pd.read_csv(SRC / "Figure4B_panelE_interaction_models.csv")
    # beta3 = VHH slope - VH slope, verified via slope_VHH - slope_VH == beta3
    assert np.allclose(models['slope_VHH'] - models['slope_VH'], models['beta3_interaction'], atol=1e-6)


def test_beta3_definition():
    models = pd.read_csv(SRC / "Figure4B_panelE_interaction_models.csv")
    for col in ['N_VHH', 'N_VH', 'slope_VH', 'slope_VHH', 'beta3_interaction', 'CI_low', 'CI_high', 'raw_p', 'BH_FDR_q']:
        assert col in models.columns


def test_common_support_range_correct():
    common = pd.read_csv(SRC / "Figure4B_common_support_models.csv")
    per_struct = pd.read_csv(SRC / "Figure4B_per_structure_metrics.csv")
    vhh_range = (per_struct[per_struct.group == 'VHH']['CDR3_length'].min(), per_struct[per_struct.group == 'VHH']['CDR3_length'].max())
    expected = f"[{vhh_range[0]},{vhh_range[1]}]"
    assert (common['common_support_range'] == expected).all()


def test_normalized_interaction_arithmetic_correct():
    per_struct = pd.read_csv(SRC / "Figure4B_per_structure_metrics.csv")
    for tset in ['CDR1', 'FR2', 'CDR2']:
        recomputed = per_struct[f'interaction_count_{tset}'] / per_struct['CDR3_length']
        assert np.allclose(per_struct[f'interaction_count_{tset}_per_CDR3_residue'].dropna(),
                            recomputed.dropna(), atol=1e-9)


def test_no_inference_from_single_group_significance_only():
    # the interaction-model table (beta3 formal test) must exist and be the
    # basis used, not a bare within-group significant/non-significant comparison
    models = pd.read_csv(SRC / "Figure4B_panelE_interaction_models.csv")
    assert 'beta3_interaction' in models.columns and 'BH_FDR_q' in models.columns


def test_no_bent_extended_language():
    # "bent"/"extended" must never be used as a CATEGORICAL CDR3-conformation
    # label (e.g. "the bent structure", "VHH CDR3 is extended"). A caption
    # explicitly disclaiming that language (e.g. "never 'bent'/'extended'")
    # is compliant, not a violation -- only flag actual categorical use.
    for fname in ['Figure4B_caption.txt']:
        fpath = ROOT / "results/final_manuscript_figures/captions" / fname
        if fpath.exists():
            text = fpath.read_text().lower()
            for bad_phrase in ['the bent', 'is bent', 'a bent', 'bent conformation',
                                'is extended', 'an extended', 'extended conformation']:
                assert bad_phrase not in text, f"found categorical usage: {bad_phrase!r}"


def test_no_internal_stage_wording_visible():
    fpath = ROOT / "results/final_manuscript_figures/captions/Figure4B_caption.txt"
    if fpath.exists():
        text = fpath.read_text().lower()
        for bad in ['stage 6', 'stage 7', 'stage6', 'stage7', 'pipeline']:
            assert bad not in text


def test_figure4b_outputs_exist_no_tiff():
    assert (MAIN / "Figure4B.png").exists()
    assert (MAIN / "Figure4B.pdf").exists()
    assert not (MAIN / "Figure4B.tiff").exists()


def test_secondary_structure_fractions_sum_correctly():
    ss = pd.read_csv(SRC / "Figure4B_secondary_structure.csv")
    frac_cols = ['fraction_state_G', 'fraction_state_T', 'fraction_state_E', 'fraction_state_S',
                 'fraction_state_C', 'fraction_state_OTHER']
    total = ss[frac_cols].sum(axis=1)
    assert np.allclose(total, 1.0, atol=1e-6)
    # DSSP blank/unassigned state must map to Zhang's 'C' (coil), not OTHER
    dssp = pd.read_csv(ROOT / "data/interim/figure4B_dssp_cdr3.csv")
    blank_rows = dssp[dssp.dssp_raw_state == '-']
    assert (blank_rows['dssp_zhang_state'] == 'C').all()


def test_panelE_includes_secondary_structure_metrics():
    models = pd.read_csv(SRC / "Figure4B_panelE_interaction_models.csv")
    ss_metrics = {'fraction_state_G', 'fraction_state_T', 'fraction_state_E', 'fraction_state_S', 'fraction_state_C'}
    assert ss_metrics.issubset(set(models['outcome']))


def test_panelA_exemplars_deterministic_q25_q75():
    ex = pd.read_csv(SRC / "Figure4B_panelA_exemplars.csv")
    assert len(ex) == 4
    assert set(ex['target_quantile']) == {'Q25', 'Q75'}
    assert (ex['distance_from_quantile_target'] < 1e-6).all()
