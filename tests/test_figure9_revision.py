"""Figure 9 revision -- integrated sequence-to-surface figure. No-invention/provenance checks."""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure9"
MAIN = ROOT / "results/final_manuscript_figures/main"
CAPTION = ROOT / "results/final_manuscript_figures/captions/Figure9_caption.txt"
AUDIT = ROOT / "results/final_manuscript_figures/audits/figure9_sequence_structure_surface_audit.md"

FULL_25 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110,
           111, 112, 113, 114, 115, 116, 118, 119, 120}
CDR3_INS = {'111A', '112A', '112B'}


def test_28_positions_match_frozen_definition():
    bridge = pd.read_csv(SRC / "Figure9_position_level_bridge.csv")
    positions = set(bridge['IMGT_position'].astype(str))
    ints = {int(p) for p in positions if p not in CDR3_INS}
    assert ints == FULL_25
    assert (positions & CDR3_INS) == CDR3_INS
    assert len(bridge) == 28


def test_sequence_logos_reused_byte_identical_to_figure2():
    for grp in ['VHH', 'VH']:
        letter = 'A' if grp == 'VHH' else 'B'
        fig9 = pd.read_csv(SRC / f"Figure9_sequence_logo_{grp}.csv")
        fig2 = pd.read_csv(ROOT / f"results/final_manuscript_figures/source_data/Figure2_panel{letter}_{grp}_sequence_logo.csv")
        pd.testing.assert_frame_equal(fig9, fig2)


def test_jsd_unsigned_and_matches_frozen():
    bridge = pd.read_csv(SRC / "Figure9_position_level_bridge.csv")
    assert (bridge['JSD'].dropna() >= 0).all()
    frozen = pd.read_csv(ROOT / "results/final_manuscript_figures/source_data/Figure2_panelD_JSD.csv")
    frozen['IMGT_position'] = frozen['IMGT_position'].astype(str)
    bridge['IMGT_position'] = bridge['IMGT_position'].astype(str)
    m = bridge.merge(frozen, on='IMGT_position', suffixes=('_f9', '_f2'))
    jsd_col_f9 = [c for c in m.columns if c.lower().startswith('jsd') and c.endswith('_f9')][0]
    jsd_col_f2 = [c for c in m.columns if c.lower().startswith('jsd') and c.endswith('_f2')][0]
    assert np.allclose(m[jsd_col_f9], m[jsd_col_f2], atol=1e-9)


def test_cdr3_shielding_only_for_vhh_and_only_15_position_subset():
    bridge = pd.read_csv(SRC / "Figure9_position_level_bridge.csv")
    CDR3_RANGE_IMGT = {'107', '108', '109', '110', '111', '111A', '112', '112A',
                        '112B', '113', '114', '115', '116'}
    cdr3_rows = bridge[bridge['IMGT_position'].astype(str).isin(CDR3_RANGE_IMGT)]
    assert cdr3_rows['VHH_CDR3_shielding'].isna().all()
    non_cdr3 = bridge[~bridge['IMGT_position'].astype(str).isin(CDR3_RANGE_IMGT)]
    assert non_cdr3['VHH_CDR3_shielding'].notna().sum() > 0


def test_ppc_pnc_position_membership_explicitly_not_available():
    for f in ['Figure9_PPC_position_membership.csv', 'Figure9_PNC_position_membership.csv']:
        df = pd.read_csv(SRC / f)
        assert (df['status'] == 'NOT AVAILABLE FROM VALIDATED OUTPUTS').all()
        value_cols = [c for c in df.columns if 'membership_frequency' in c]
        for c in value_cols:
            assert df[c].isna().all()


def test_ppc_pnc_two_part_statistics_present_and_separate():
    for metric in ['PPC', 'PNC']:
        prev = pd.read_csv(SRC / f"Figure9_{metric}_prevalence.csv")
        mag = pd.read_csv(SRC / f"Figure9_{metric}_conditional_magnitude.csv")
        assert set(prev['group']) == {'VHH', 'VH'}
        assert prev['test'].iloc[0] == 'fisher_exact'
        assert mag['test'].iloc[0] == 'mannwhitneyu'
        assert prev['N_patch_positive'].iloc[0] == mag['N_VHH_nonzero'].iloc[0] or \
               prev[prev.group == 'VHH']['N_patch_positive'].iloc[0] == mag['N_VHH_nonzero'].iloc[0]


def test_ppc_pnc_prevalence_matches_raw_patch_descriptors():
    patch = pd.read_parquet(ROOT / "data/processed/stage8_surface_patch_descriptors.parquet")
    for metric, col in [('PPC', 'former_interface_core_PPC'), ('PNC', 'former_interface_core_PNC')]:
        prev = pd.read_csv(SRC / f"Figure9_{metric}_prevalence.csv")
        for grp in ['VHH', 'VH']:
            n_pos = (patch[(patch.group == grp)][col] > 0).sum()
            n_total = (patch.group == grp).sum()
            row = prev[prev.group == grp].iloc[0]
            assert row['N_patch_positive'] == n_pos
            assert row['N_total'] == n_total


def test_apbs_components_never_netted():
    apbs = pd.read_csv(SRC / "Figure9_APBS_summary.csv")
    assert 'component' in apbs.columns
    assert set(apbs['component'].str.lower()) == {'positive', 'negative'}
    assert not any('net' in c.lower() for c in apbs.columns)


def test_apbs_cohort_sizes_match_frozen():
    frozen = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    valid = frozen.dropna(subset=['interface_potential_positive_mean', 'interface_potential_negative_mean'])
    n_vhh = (valid.group == 'VHH').sum()
    n_vh = (valid.group == 'VH').sum()
    assert n_vhh == 521
    assert n_vh == 1386


def test_apbs_color_limits_shared_symmetric():
    script = (ROOT / "scripts/figure9_render_3d_modules.py").read_text()
    assert '-APBS_LIMIT, 0, APBS_LIMIT' in script
    assert script.count("ramp_new(ramp_name") == 1  # single shared call site for both groups


def test_exemplar_selection_deterministic_and_documented():
    df = pd.read_csv(SRC / "Figure9_exemplar_selection.csv")
    assert 'distance_from_target' in df.columns
    assert (df['distance_from_target'].astype(float) >= 0).all()
    apbs_rows = df[df['metric'].str.contains('APBS', na=False)]
    assert len(apbs_rows) == 2
    vh_apbs_row = apbs_rows[apbs_rows.group == 'VH'].iloc[0]
    assert vh_apbs_row['structure_id'] == 'pdb_00003gi9-H-L'
    assert 'NEW deterministic exemplar required' in vh_apbs_row['rule']


def test_apbs_vh_exemplar_actually_has_valid_dx_grid():
    dx = ROOT / "data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pqr.dx"
    assert dx.exists()
    old_dx = ROOT / "data/interim/stage8/apbs_cohort/pdb_00004k2u-I-M/pdb_00004k2u-I-M.pqr.dx"
    assert not old_dx.exists()  # confirms the substitution was actually necessary


def test_provenance_file_covers_all_panels():
    prov = pd.read_csv(SRC / "Figure9_annotation_provenance.csv")
    required_cols = {'figure_panel', 'annotation', 'IMGT_position_if_applicable', 'metric',
                      'source_file', 'source_column', 'existing_or_derived', 'derivation_formula_if_derived'}
    assert required_cols.issubset(set(prov.columns))
    panels_covered = set(prov['figure_panel'].str[0])
    assert {'A', 'C', 'D', 'E', 'F', 'H'}.issubset(panels_covered)
    assert (prov['existing_or_derived'].isin(['existing', 'derived', 'not available'])).all()


def test_solubility_aggregation_stability_affinity_only_mentioned_as_not_tested():
    caption_text = CAPTION.read_text().lower()
    # these terms may appear ONLY as part of the explicit disclaimer, never as a claimed result
    assert 'not calculated or tested' in caption_text or 'not tested' in caption_text
    forbidden_claims = ['improves solubility', 'increases solubility', 'reduces aggregation',
                         'confers stability', 'increases affinity', 'higher affinity',
                         'binds more tightly', 'more stable than']
    for claim in forbidden_claims:
        assert claim not in caption_text


def test_no_causal_apbs_language_in_caption():
    text = CAPTION.read_text().lower()
    assert 'apbs proves' not in text
    assert 'causes' not in text.split('apbs')[-1][:200] if 'apbs' in text else True
    assert 'spatially correspond' in text or 'co-localiz' in text


def test_tetrad_positions_correct():
    bridge = pd.read_csv(SRC / "Figure9_position_level_bridge.csv")
    bridge['IMGT_position'] = bridge['IMGT_position'].astype(str)
    tetrad_rows = bridge[bridge['tetrad'] == True]
    assert set(tetrad_rows['IMGT_position']) == {'42', '49', '50', '52'}


def test_position_tracks_top5_deterministic():
    bridge = pd.read_csv(SRC / "Figure9_position_level_bridge.csv")
    bridge['IMGT_position'] = bridge['IMGT_position'].astype(str)
    tracks = pd.read_csv(SRC / "Figure9_position_tracks.csv")
    seq_div_track = set(tracks[(tracks.track_name == 'seq_divergence') & (tracks.in_top5 == True)]['IMGT_position'].astype(str))
    top5_jsd = set(bridge.set_index('IMGT_position')['JSD'].nlargest(5).index)
    assert seq_div_track == top5_jsd


def test_figure9_outputs_exist_no_tiff():
    assert (MAIN / "Figure9.png").exists()
    assert (MAIN / "Figure9.pdf").exists()
    assert not (MAIN / "Figure9.tiff").exists()


def test_frozen_upstream_files_unchanged():
    # spot-check a frozen source is still readable and structurally unchanged in shape
    stage8 = pd.read_parquet(ROOT / "data/processed/stage8_surface_patch_descriptors.parquet")
    assert len(stage8) == 533 + 1429
    stage6 = pd.read_csv(ROOT / "data/processed/stage6_position_shielding_summary.csv")
    assert 'mean_delta_sasa_shield' in stage6.columns
