"""View-only helper: interactively inspect ONE Figure4B CDR3-geometry
exemplar at a time (env vhh-interface-viz). NOT a modification of
scripts/figure4B_render_panelA.py (the original 4-panel version -- this is the
one actually described in Figure4B_caption.txt: "CDR1 (gold), FR2 (teal),
CDR2 (green), CDR3 (red), IMGT42 (magenta sphere)", not the alternate
Zhang-color-matched figure4B_render_panelA_v2.py) -- duplicates its exact
load/color logic (unchanged) for a single chosen structure, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select with env var FIG4B_WHICH (default: vhh_short):
  vhh_short / vhh_long -- VHH, Q25/Q75 CDR3-length exemplar
  vh_short  / vh_long  -- conventional VH, Q25/Q75 CDR3-length exemplar (VL translucent)
"""
from pymol import cmd
from pathlib import Path
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure4b_structures"

CDR3_RANGE = list(range(105, 118))
CDR1_RANGE = list(range(27, 39))
FR2_RANGE = list(range(39, 56))
CDR2_RANGE = list(range(56, 66))

COL = {'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
       'cdr3_red': [0.85, 0.15, 0.15], 'cdr1_gold': [0.72, 0.53, 0.04], 'fr2_teal': [0.0, 0.55, 0.55],
       'cdr2_green': [0.35, 0.65, 0.15], 'pos42_magenta': [0.85, 0.1, 0.65]}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)


def color_regions(obj, chain, main_color):
    sel = f'{obj} and chain {chain}'
    cmd.color(main_color, sel)
    cmd.color('cdr1_gold', f'{sel} and resi {"+".join(map(str, CDR1_RANGE))}')
    cmd.color('fr2_teal', f'{sel} and resi {"+".join(map(str, FR2_RANGE))}')
    cmd.color('cdr2_green', f'{sel} and resi {"+".join(map(str, CDR2_RANGE))}')
    cmd.color('cdr3_red', f'{sel} and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.show('spheres', f'{sel} and resi 42 and name CA')
    cmd.color('pos42_magenta', f'{sel} and resi 42 and name CA')
    cmd.set('sphere_scale', 0.7, obj)


which = os.environ.get('FIG4B_WHICH', 'vhh_short')
chain, color = (('A', 'vhh_blue') if which.startswith('vhh') else ('H', 'vh_purple'))

cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('antialias', 2)

cmd.load(str(STRUCT_DIR / f'{which}.pdb'), which)
cmd.hide('everything', which)
cmd.show('cartoon', f'{which} and chain {chain}')
color_regions(which, chain, color)
if which.startswith('vh_'):
    cmd.show('cartoon', f'{which} and chain L')
    cmd.color('vl_pink', f'{which} and chain L')
    cmd.set('cartoon_transparency', 0.7, f'{which} and chain L')

cmd.orient(f'{which} and chain {chain}')
cmd.zoom(f'{which} and chain {chain}', buffer=3)
print(f'Figure4B {which} -- gold=CDR1, teal=FR2, green=CDR2, red=CDR3, magenta sphere=IMGT42.')
