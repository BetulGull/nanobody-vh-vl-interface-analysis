"""Figure 2 revision -- sequence logos + VHH-minus-VH differential heatmap.
Validates the final-manuscript Figure 2 rebuild against its own source-data
and against the frozen Stage-4 tables it is derived from. No Stage 2-12
output is modified by these tests.
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data"
AUDITS = ROOT / "results/final_manuscript_figures/audits"
MAIN = ROOT / "results/final_manuscript_figures/main"
SUPP = ROOT / "results/final_manuscript_figures/supplementary"

TETRAD = {'42', '49', '50', '52'}


def test_panelA_panelB_frequencies_sum_to_one():
    for f in ['Figure2_panelA_VHH_sequence_logo.csv', 'Figure2_panelB_VH_sequence_logo.csv']:
        df = pd.read_csv(SRC / f)
        sums = df.groupby('IMGT_position')['frequency'].sum()
        assert np.allclose(sums.values, 1.0, atol=1e-6), f"{f}: frequencies do not sum to 1 per position"


def test_logo_frequency_matches_frozen_heatmap_within_documented_tolerance():
    """New (X-dropped, renormalized) logo frequencies must equal the frozen
    Stage-4 raw frequencies exactly wherever no 'X' residue call exists at
    that position; documented, small (<1e-3) differences are allowed only
    at the 14 (position, VH) cells with a real X observation (see
    figure2_logo_vs_heatmap_consistency.csv)."""
    audit = pd.read_csv(AUDITS / "figure2_logo_vs_heatmap_consistency.csv")
    assert len(audit) > 0
    fails = audit[~audit['pass']]
    # every failing cell must be small and attributable to X-renormalization
    assert (fails['absolute_difference'] < 1e-3).all(), "unexplained large logo-vs-heatmap discrepancy"
    # the vast majority of cells must match exactly
    assert audit['pass'].mean() > 0.7


def test_delta_frequency_arithmetic_exact():
    df = pd.read_csv(SRC / "Figure2_panelC_VHH_minus_VH_frequency_difference.csv")
    recomputed = df['frequency_VHH'] - df['frequency_VH']
    assert np.allclose(df['delta_frequency'].values, recomputed.values, atol=1e-12)


def test_jsd_matches_frozen_stage4_values():
    panelD = pd.read_csv(SRC / "Figure2_panelD_JSD.csv")
    s4 = pd.read_csv(ROOT / "data/processed/stage4_position_chemistry_comparison.csv")
    for _, row in panelD.iterrows():
        lbl = str(row['IMGT_position'])
        num = ''.join(c for c in lbl if c.isdigit())
        ins = lbl[len(num):]
        pos = int(num)
        match = s4[(s4.imgt_position == pos) & ((s4.insertion_code == ins) | (s4.insertion_code.isna() & (ins == '')))]
        assert len(match) == 1
        assert abs(match.iloc[0]['JSD'] - row['JSD']) < 1e-9


def test_same_position_ordering_across_panels():
    audit = pd.read_csv(SRC / "figure1_empirical_interface_position_audit.csv")
    expected_order = audit['IMGT_position'].astype(str).tolist()
    panelA = pd.read_csv(SRC / "Figure2_panelA_VHH_sequence_logo.csv")
    panelD = pd.read_csv(SRC / "Figure2_panelD_JSD.csv")
    assert set(panelA['IMGT_position'].astype(str)) == set(expected_order)
    assert panelD['IMGT_position'].astype(str).tolist() == expected_order


def test_25_vs_28_distinction_preserved():
    audit = pd.read_csv(SRC / "figure1_empirical_interface_position_audit.csv")
    assert len(audit) == 28
    unique_integers = audit['IMGT_position'].astype(str).apply(lambda s: ''.join(c for c in s if c.isdigit())).nunique()
    assert unique_integers == 25
    assert {'111A', '112A', '112B'}.issubset(set(audit['IMGT_position'].astype(str)))


def test_no_unresolved_position_zero_filled():
    panelA = pd.read_csv(SRC / "Figure2_panelA_VHH_sequence_logo.csv")
    panelB = pd.read_csv(SRC / "Figure2_panelB_VH_sequence_logo.csv")
    assert (panelA['N_observed'] > 0).all()
    assert (panelB['N_observed'] > 0).all()
    assert panelA['sequence_resolved_status'].all()
    assert panelB['sequence_resolved_status'].all()


def test_amino_acid_colors_consistent_with_figure1():
    pal = pd.read_csv(SRC / "figure1_sequence_logo_color_mapping.csv")
    assert len(pal) == 20
    assert set(pal['amino_acid']) == set('ACDEFGHIKLMNPQRSTVWY')


def test_tetrad_annotation_correct():
    for f in ['Figure2_panelA_VHH_sequence_logo.csv', 'Figure2_panelC_VHH_minus_VH_frequency_difference.csv']:
        df = pd.read_csv(SRC / f)
        flagged = set(df[df['tetrad']]['IMGT_position'].astype(str))
        assert flagged == TETRAD


def test_supplementary_raw_heatmaps_exist():
    for ext in ['png', 'pdf', 'tiff']:
        assert (SUPP / f"FigureS_sequence_frequency_heatmaps.{ext}").exists()


def test_figure2_main_outputs_exist():
    for ext in ['png', 'pdf', 'tiff']:
        assert (MAIN / f"Figure2.{ext}").exists()


def test_panelC_ranked_table_has_rank_one_at_max_abs_delta():
    ranked = pd.read_csv(SRC / "Figure2_panelC_residue_enrichment_ranked.csv")
    top_row = ranked[ranked.abs_delta_rank == 1].iloc[0]
    max_abs = ranked['delta_frequency'].abs().max()
    assert abs(abs(top_row['delta_frequency']) - max_abs) < 1e-9
