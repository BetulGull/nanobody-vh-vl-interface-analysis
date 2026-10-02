"""View-only helper: interactively inspect ONE Figure7 complex at a time (env
vhh-interface-viz). NOT a modification of scripts/figure7_render_revised.py --
duplicates its exact load/color/patch logic (unchanged) for a single chosen
complex, skips the final cmd.ray()/cmd.png() save step, and leaves the PyMOL
GUI open so it can be rotated/zoomed by hand.

Select with env var FIG7_WHICH (default: vhh):
  vhh -- VHH + antigen complex (panel B)
  vh  -- conventional VH + VL + antigen complex (panel C)
"""
from pymol import cmd
from pathlib import Path
import csv
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure7_structures"
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure7"

COL = {
    'vhh_blue': [0.19, 0.38, 0.80], 'vh_purple': [0.49, 0.31, 0.66], 'vl_pink': [0.95, 0.65, 0.77],
    'antigen_gray': [0.78, 0.78, 0.78], 'epitope_darkgray': [0.44, 0.44, 0.44],
    'fi_yellow': [1.0, 0.83, 0.0], 'fi_red': [0.82, 0.08, 0.08], 'other_orange': [1.0, 0.60, 0.0],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)

MATCHED_EPITOPE_RESI = [41, 42, 43, 44, 45, 47, 66, 67, 69, 71, 72]


def apply_residue_classes(obj, chain, csv_path):
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    for cls, color in [('former_interface_nonengaged', 'fi_yellow'),
                        ('former_interface_engaged', 'fi_red'),
                        ('other_antigen_contact', 'other_orange')]:
        positions = [r['IMGT_position'] for r in rows if r['render_class'] == cls]
        if positions:
            sel = '+'.join(str(p) for p in positions)
            patch_sel = f'{obj} and chain {chain} and resi {sel}'
            cmd.show('surface', patch_sel)
            cmd.color(color, patch_sel)


which = os.environ.get('FIG7_WHICH', 'vhh')
cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('antialias', 2)
cmd.set('cartoon_transparency', 0.0)
cmd.set('surface_quality', 1)

epi_sel = '+'.join(str(r) for r in MATCHED_EPITOPE_RESI)

if which == 'vhh':
    cmd.load(str(STRUCT_DIR / 'vhh_complex_antigen_aligned.pdb'), 'vhh')
    cmd.hide('everything', 'vhh')
    cmd.show('cartoon', 'vhh and chain A')
    cmd.color('vhh_blue', 'vhh and chain A')
    cmd.show('surface', 'vhh and chain G')
    cmd.color('antigen_gray', 'vhh and chain G')
    cmd.color('epitope_darkgray', f'vhh and chain G and resi {epi_sel}')
    apply_residue_classes('vhh', 'A', SRC / 'Figure7_structure_B_residue_classes.csv')
    cmd.orient('vhh')
    cmd.zoom('vhh', buffer=5)
    print('Figure7 VHH complex (panel B) -- blue cartoon VHH, gray surface antigen, '
          'dark-gray matched epitope, yellow/red/orange classified residue patches.')
else:
    cmd.load(str(STRUCT_DIR / 'vh_complex_antigen_aligned.pdb'), 'vh')
    cmd.hide('everything', 'vh')
    cmd.show('cartoon', 'vh and chain H')
    cmd.color('vh_purple', 'vh and chain H')
    cmd.show('cartoon', 'vh and chain L')
    cmd.color('vl_pink', 'vh and chain L')
    cmd.show('surface', 'vh and chain G')
    cmd.color('antigen_gray', 'vh and chain G')
    cmd.color('epitope_darkgray', f'vh and chain G and resi {epi_sel}')
    apply_residue_classes('vh', 'H', SRC / 'Figure7_structure_C_residue_classes.csv')
    cmd.orient('vh')
    cmd.zoom('vh', buffer=5)
    print('Figure7 VH+VL complex (panel C) -- purple cartoon VH, pink cartoon VL, '
          'gray surface antigen, dark-gray matched epitope, yellow/red/orange classified residue patches.')
