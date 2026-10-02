"""View-only helper: interactively inspect ONE Figure6 CDR3-packing exemplar
at a time (env vhh-interface-viz). NOT a modification of
scripts/final_manuscript_render_structures.py::figure6() -- duplicates its
exact load/color logic (unchanged) for a single chosen structure, skips the
final cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can
be rotated/zoomed by hand. VHH-only figure (no conventional-VH counterpart).

Select with env var FIG6_WHICH (default: short):
  short -- VHH, nearest-to-Q25 CDR3 length
  long  -- VHH, nearest-to-Q75 CDR3 length
Both: blue body, teal=CDR3, orange=14-position former-interface packing target.
"""
from pymol import cmd
from pathlib import Path
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/final_manuscript_structures"

FULL_28 = [40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120]
CDR3_RANGE = list(range(105, 118))

COL = {
    'vhh_blue': [0.20, 0.40, 0.80], 'cdr3_teal': [0.0, 0.60, 0.60],
    'interface_orange': [1.0, 0.60, 0.0],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)

which = os.environ.get('FIG6_WHICH', 'short')
fname = 'fig6_vhh_short.pdb' if which == 'short' else 'fig6_vhh_long_aligned.pdb'
obj = f'vhh_{which}'

cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('cartoon_fancy_helices', 1)
cmd.set('antialias', 2)

cmd.load(str(STRUCT_DIR / fname), obj)
cmd.hide('everything', obj)
cmd.show('cartoon', obj)
cmd.color('vhh_blue', obj)
cmd.color('cdr3_teal', f'{obj} and resi {"+".join(map(str, CDR3_RANGE))}')
cmd.color('interface_orange', f'{obj} and resi {"+".join(map(str, FULL_28))} and not resi {"+".join(map(str, CDR3_RANGE))}')

cmd.orient(obj)
cmd.zoom(obj, buffer=3)
label = 'Q25 (short)' if which == 'short' else 'Q75 (long)'
print(f'Figure6 VHH {label} CDR3-length exemplar -- blue body, teal=CDR3, orange=14-position packing target.')
