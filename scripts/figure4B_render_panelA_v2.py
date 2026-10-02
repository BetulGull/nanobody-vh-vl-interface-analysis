"""Figure 4B Panel A -- CDR3 geometry, TWO structures (VHH short vs VHH long CDR3),
matching Zhang Figure 4A's own single-group illustration exactly (their panel is not
species-split either). PyMOL script (env vhh-interface-viz). Reuses the existing
deterministic Q25/Q75 VHH exemplars frozen in Figure4B_panelA_exemplars.csv and the
already-aligned structure files in data/interim/figure4b_structures/ -- no new
structure selection, no recomputation.

Color mapping matches the Zhang reference caption exactly:
CDR3=red, CDR1=orange, FR2=cyan, CDR2=magenta, IMGT42=pink marker.
"""
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure4b_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"

CDR3_RANGE = list(range(105, 118))
CDR1_RANGE = list(range(27, 39))
FR2_RANGE = list(range(39, 56))
CDR2_RANGE = list(range(56, 66))

COL = {
    'vhh_blue': [0.20, 0.40, 0.80],
    'cdr3_red': [0.85, 0.15, 0.15],
    'cdr1_orange': [1.0, 0.60, 0.0],
    'fr2_cyan': [0.0, 0.75, 0.85],
    'cdr2_magenta': [0.80, 0.15, 0.65],
    'pos42_pink': [0.95, 0.55, 0.70],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)


def base_settings():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('ambient', 0.4)
    cmd.set('specular', 0.15)
    cmd.set('antialias', 2)


def ray_png(path, width=1200, height=1200, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


def color_regions(obj, chain):
    sel = f'{obj} and chain {chain}'
    cmd.color('vhh_blue', sel)
    cmd.color('cdr1_orange', f'{sel} and resi {"+".join(map(str, CDR1_RANGE))}')
    cmd.color('fr2_cyan', f'{sel} and resi {"+".join(map(str, FR2_RANGE))}')
    cmd.color('cdr2_magenta', f'{sel} and resi {"+".join(map(str, CDR2_RANGE))}')
    cmd.color('cdr3_red', f'{sel} and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.show('spheres', f'{sel} and resi 42 and name CA')
    cmd.color('pos42_pink', f'{sel} and resi 42 and name CA')
    cmd.set('sphere_scale', 0.7, obj)


cmd.delete('all')
base_settings()
for name in ['vhh_short', 'vhh_long']:
    cmd.load(str(STRUCT_DIR / f'{name}.pdb'), name)
    cmd.hide('everything', name)
    cmd.show('cartoon', f'{name} and chain A')
    color_regions(name, 'A')
    cmd.orient(f'{name} and chain A')
    cmd.zoom(f'{name} and chain A', buffer=3)
    ray_png(OUT_DIR / f'Figure4B_panelA_v2_{name}.png')

print('Figure 4B Panel A (2 structures, Zhang color mapping) rendered')
