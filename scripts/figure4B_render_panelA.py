"""Figure 4B Panel A -- CDR3 geometry, 2x2 (VHH/VH x short/long CDR3).
PyMOL script (env vhh-interface-viz). Deterministic Q25/Q75 exemplars
frozen in Figure4B_panelA_exemplars.csv, common-frame aligned."""
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure4b_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"

CDR3_RANGE = list(range(105, 118))
CDR1_RANGE = list(range(27, 39))
FR2_RANGE = list(range(39, 56))
CDR2_RANGE = list(range(56, 66))

COL = {'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
       'cdr3_red': [0.85, 0.15, 0.15], 'cdr1_gold': [0.72, 0.53, 0.04], 'fr2_teal': [0.0, 0.55, 0.55],
       'cdr2_green': [0.35, 0.65, 0.15], 'pos42_magenta': [0.85, 0.1, 0.65]}
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


cmd.delete('all')
base_settings()
structs = {'vhh_short': ('A', 'vhh_blue'), 'vhh_long': ('A', 'vhh_blue'),
           'vh_short': ('H', 'vh_purple'), 'vh_long': ('H', 'vh_purple')}
for name, (chain, color) in structs.items():
    cmd.load(str(STRUCT_DIR / f'{name}.pdb'), name)
    cmd.hide('everything', name)
    cmd.show('cartoon', f'{name} and chain {chain}')
    color_regions(name, chain, color)
    if 'vh_' in name and 'vhh' not in name:
        cmd.show('cartoon', f'{name} and chain L')
        cmd.color('vl_pink', f'{name} and chain L')
        cmd.set('cartoon_transparency', 0.7, f'{name} and chain L')

cmd.orient('vhh_short')
cmd.zoom('vhh_short', buffer=3)
view = cmd.get_view()

for name in structs:
    for other in structs:
        cmd.disable(other)
    cmd.enable(name)
    cmd.set_view(view)
    ray_png(OUT_DIR / f'Figure4B_panelA_{name}.png')

print('Figure 4B Panel A rendered (4 panels)')
