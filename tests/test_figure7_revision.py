"""Figure 7 revision -- explicit 28-position former-interface engagement coloring."""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure7"
MAIN = ROOT / "results/final_manuscript_figures/main"
CAPTION = ROOT / "results/final_manuscript_figures/captions/Figure7_caption.txt"

FULL_25 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110,
           111, 112, 113, 114, 115, 116, 118, 119, 120}


def test_28_position_set_matches_frozen_definition():
    for f in ['Figure7_structure_B_residue_classes.csv', 'Figure7_structure_C_residue_classes.csv']:
        df = pd.read_csv(SRC / f)
        fi_positions = set(df[df.former_interface_28]['IMGT_position'])
        assert fi_positions == FULL_25


def test_engaged_subset_uses_same_contact_definition_as_panel_a():
    eng = pd.read_csv(SRC / "Figure7_cohort_engagement.csv")
    b_frozen = eng[eng.INSTANCE == 'pdb_00008k4q-B'].iloc[0]
    c_frozen = eng[eng.INSTANCE == 'pdb_00006wgl-A-B'].iloc[0]

    df_b = pd.read_csv(SRC / "Figure7_structure_B_residue_classes.csv")
    df_c = pd.read_csv(SRC / "Figure7_structure_C_residue_classes.csv")
    n_engaged_b = (df_b.render_class == 'former_interface_engaged').sum()
    n_engaged_c = (df_c.render_class == 'former_interface_engaged').sum()

    assert n_engaged_b == b_frozen['full_n_contacting']
    assert n_engaged_c == c_frozen['full_n_contacting']


def test_classification_mutually_exclusive():
    for f in ['Figure7_structure_B_residue_classes.csv', 'Figure7_structure_C_residue_classes.csv']:
        df = pd.read_csv(SRC / f)
        assert df['render_class'].isin(
            ['former_interface_engaged', 'former_interface_nonengaged', 'other_antigen_contact']
        ).all()
        dup = df.groupby('IMGT_position').size()
        assert (dup == 1).all()


def test_no_former_interface_residue_colored_orange():
    for f in ['Figure7_structure_B_residue_classes.csv', 'Figure7_structure_C_residue_classes.csv']:
        df = pd.read_csv(SRC / f)
        fi_rows = df[df.former_interface_28]
        assert (fi_rows['render_class'] != 'other_antigen_contact').all()


def test_no_non_former_interface_residue_colored_red():
    for f in ['Figure7_structure_B_residue_classes.csv', 'Figure7_structure_C_residue_classes.csv']:
        df = pd.read_csv(SRC / f)
        non_fi = df[~df.former_interface_28]
        assert (non_fi['render_class'] != 'former_interface_engaged').all()


def test_yellow_equals_F_minus_A_red_equals_F_intersect_A():
    for f in ['Figure7_structure_B_residue_classes.csv', 'Figure7_structure_C_residue_classes.csv']:
        df = pd.read_csv(SRC / f)
        fi = df[df.former_interface_28]
        yellow = set(fi[fi.render_class == 'former_interface_nonengaged']['IMGT_position'])
        red = set(fi[fi.render_class == 'former_interface_engaged']['IMGT_position'])
        engaged_positions = set(fi[fi.antigen_contact]['IMGT_position'])
        assert red == engaged_positions
        assert yellow == (FULL_25 - engaged_positions)


def test_color_mapping_antigen_is_gray_not_green():
    cmap = pd.read_csv(SRC / "Figure7_color_mapping.csv")
    antigen_row = cmap[cmap.element == 'Antigen body'].iloc[0]
    r, g, b = [int(antigen_row['hex'][i:i+2], 16) for i in (1, 3, 5)]
    # gray = R≈G≈B, not green-dominant
    assert abs(r - g) < 15 and abs(g - b) < 15
    assert 'CHANGED from green' in antigen_row['notes']


def test_boxplot_median_not_orange():
    cmap = pd.read_csv(SRC / "Figure7_color_mapping.csv")
    median_row = cmap[cmap.element.str.contains('median', case=False)].iloc[0]
    assert median_row['hex'].upper() != '#FFA500'
    assert 'CHANGED from orange' in median_row['notes']


def test_vhh_box_color_matches_vhh_structural_family():
    cmap = pd.read_csv(SRC / "Figure7_color_mapping.csv")
    vhh_struct = cmap[cmap.element == 'VHH domain'].iloc[0]['hex']
    vhh_box = cmap[cmap.element.str.contains('VHH box')].iloc[0]['hex']
    def to_rgb(h):
        return tuple(int(h[i:i+2], 16) for i in (1, 3, 5))
    s, b = to_rgb(vhh_struct), to_rgb(vhh_box)
    # same hue family: blue channel dominant in both
    assert s[2] > s[0] and b[2] > b[0]


def test_conventional_box_derives_from_vh_vl_family():
    cmap = pd.read_csv(SRC / "Figure7_color_mapping.csv")
    box_row = cmap[cmap.element.str.contains('conventional VH box', case=False)].iloc[0]
    assert 'lavender' in box_row['notes'].lower()
    assert 'VH purple' in box_row['notes'] or 'purple' in box_row['notes'].lower()


def test_matched_epitope_only_if_validated():
    mep = pd.read_csv(SRC / "Figure7_matched_pair_summary.csv")
    assert (mep['status'] == 'VALID_MATCHED_EPITOPE').all()


def test_underlying_engagement_values_unchanged():
    eng = pd.read_csv(SRC / "Figure7_cohort_engagement.csv")
    frozen = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    frozen = frozen[['INSTANCE', 'full_interface_engagement_fraction']]
    m = eng.merge(frozen, on='INSTANCE', suffixes=('_new', '_frozen'))
    valid = m.dropna(subset=['full_interface_engagement_fraction_new', 'full_interface_engagement_fraction_frozen'])
    assert np.allclose(valid['full_interface_engagement_fraction_new'],
                        valid['full_interface_engagement_fraction_frozen'], atol=1e-9)


def test_no_stage_number_in_caption_or_titles():
    text = CAPTION.read_text().lower()
    assert 'stage 3' not in text and 'stage3' not in text
    assert 'stage 9' not in text and 'stage9' not in text


def test_figure7_outputs_exist_no_tiff():
    assert (MAIN / "Figure7.png").exists()
    assert (MAIN / "Figure7.pdf").exists()
    assert not (MAIN / "Figure7.tiff").exists()


def test_vl_never_classified_as_former_interface():
    df_c = pd.read_csv(SRC / "Figure7_structure_C_residue_classes.csv")
    assert (df_c['antibody_chain'] != 'VL').all()
