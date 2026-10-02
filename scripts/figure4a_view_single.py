"""View-only helper: interactively inspect ONE Figure4A structural-context
exemplar at a time (env vhh-interface-viz). NOT a modification of
scripts/figure4A_render_structural_context.py -- duplicates its exact
load/color logic (unchanged) for a single chosen structure, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select with env var FIG4A_WHICH (default: vhh):
  vhh -- VHH exemplar pdb_00008qf4-E
  vh  -- conventional VH exemplar pdb_00009f1i-H-L (VL translucent for context)
Both: red=CDR3, orange=15-position shieldable former-interface, magenta
spheres=FR2 tetrad (42/49/50/52).
"""
from pymol import cmd
from pathlib import Path
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure4a_structures"

SHIELDABLE_15 = [40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120]
CDR3_RANGE = list(range(105, 118))
TETRAD = [42, 49, 50, 52]

COL = {'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
       'cdr3_red': [0.85, 0.15, 0.15], 'iface_orange': [1.0, 0.60, 0.0], 'tetrad_magenta': [0.85, 0.1, 0.65]}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)

which = os.environ.get('FIG4A_WHICH', 'vhh')

cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('antialias', 2)

if which == 'vhh':
    cmd.load(str(STRUCT_DIR / 'vhh_exemplar.pdb'), 'vhh')
    cmd.hide('everything', 'vhh')
    cmd.show('cartoon', 'vhh and chain A')
    cmd.color('vhh_blue', 'vhh and chain A')
    cmd.color('iface_orange', f'vhh and chain A and resi {"+".join(map(str, SHIELDABLE_15))}')
    cmd.color('cdr3_red', f'vhh and chain A and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.show('spheres', f'vhh and chain A and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.color('tetrad_magenta', f'vhh and chain A and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.set('sphere_scale', 0.6, 'vhh')
    cmd.orient('vhh')
    cmd.zoom('vhh', buffer=3)
    print('Figure4A VHH exemplar (pdb_00008qf4-E) -- blue body, red=CDR3, '
          'orange=15-position shieldable interface, magenta spheres=FR2 tetrad.')
else:
    cmd.load(str(STRUCT_DIR / 'vh_exemplar_aligned.pdb'), 'vh')
    cmd.hide('everything', 'vh')
    cmd.show('cartoon', 'vh and chain H')
    cmd.color('vh_purple', 'vh and chain H')
    cmd.color('iface_orange', f'vh and chain H and resi {"+".join(map(str, SHIELDABLE_15))}')
    cmd.color('cdr3_red', f'vh and chain H and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.show('spheres', f'vh and chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.color('tetrad_magenta', f'vh and chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.set('sphere_scale', 0.6, 'vh')
    cmd.show('cartoon', 'vh and chain L')
    cmd.color('vl_pink', 'vh and chain L')
    cmd.set('cartoon_transparency', 0.7, 'vh and chain L')
    cmd.orient('vh')
    cmd.zoom('vh', buffer=3)
    print('Figure4A conventional-VH exemplar (pdb_00009f1i-H-L) -- purple VH body, '
          'red=CDR3, orange=15-position shieldable interface, magenta spheres=FR2 tetrad, pink VL translucent.')
