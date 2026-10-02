"""View-only helper: interactively inspect Figure1's native Fv structure, at
either coloring stage (env vhh-interface-viz). NOT a modification of
scripts/final_manuscript_render_structures.py::figure1() -- duplicates its
exact load/color logic (unchanged) for a single chosen stage, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select with env var FIG1_WHICH (default: panelb):
  panela -- Panel A state: VH (purple) + VL (pink) cartoon, CDR3 highlighted teal
  panelb -- Panel B state (final): VH former-interface (orange), shieldable-15
            (yellow), FR2 tetrad (magenta spheres), VL faded/translucent
"""
from pymol import cmd
from pathlib import Path
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/final_manuscript_structures"

FULL_28 = [40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120]
CDR3_RANGE = list(range(105, 118))
SHIELDABLE_15 = [p for p in FULL_28 if p not in CDR3_RANGE]
TETRAD = [42, 49, 50, 52]

COL = {
    'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
    'antigen_green': [0.55, 0.75, 0.55], 'cdr3_teal': [0.0, 0.60, 0.60],
    'interface_orange': [1.0, 0.60, 0.0], 'shieldable_yellow': [1.0, 0.85, 0.2],
    'tetrad_magenta': [0.85, 0.1, 0.65], 'grey': [0.82, 0.82, 0.82],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)

which = os.environ.get('FIG1_WHICH', 'panelb')

cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('cartoon_fancy_helices', 1)
cmd.set('antialias', 2)

cmd.load(str(STRUCT_DIR / 'figure1_native_fv.pdb'), 'fv')
cmd.hide('everything')
cmd.show('cartoon', 'fv')
cmd.color('vh_purple', 'chain H')
cmd.color('vl_pink', 'chain L')
cmd.set('cartoon_transparency', 0.0)
cmd.orient('fv')
cmd.zoom('fv', buffer=3)
cmd.turn('y', 20)

cmd.color('cdr3_teal', f'chain H and resi {"+".join(map(str, CDR3_RANGE))}')

if which == 'panelb':
    cmd.color('vh_purple', 'chain H')
    cmd.color('interface_orange', f'chain H and resi {"+".join(map(str, FULL_28))}')
    cmd.color('shieldable_yellow', f'chain H and resi {"+".join(map(str, SHIELDABLE_15))}')
    cmd.show('spheres', f'chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.color('tetrad_magenta', f'chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.set('sphere_scale', 0.6)
    cmd.set('cartoon_transparency', 0.65, 'chain L')
    print('Figure1 Panel B state -- purple VH, orange=28-pos former interface, '
          'yellow=shieldable-15, magenta spheres=FR2 tetrad, pink VL faded.')
else:
    print('Figure1 Panel A state -- purple VH, pink VL, teal=CDR3.')
