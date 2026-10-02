"""Figure 4A Panel A -- structural context. PyMOL script (env vhh-interface-viz).
Deterministic exemplars: pdb_00008qf4-E (VHH), pdb_00009f1i-H-L (VH), frozen
in Figure4A_structure_exemplars.csv, common-framework-aligned (Figure4A_alignment_audit.csv)."""
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure4a_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"

SHIELDABLE_15 = [40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120]
CDR3_RANGE = list(range(105, 118))
TETRAD = [42, 49, 50, 52]

COL = {'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
       'cdr3_red': [0.85, 0.15, 0.15], 'iface_orange': [1.0, 0.60, 0.0], 'tetrad_magenta': [0.85, 0.1, 0.65]}
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


def ray_png(path, width=1400, height=1400, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


cmd.delete('all')
base_settings()
cmd.load(str(STRUCT_DIR / 'vhh_exemplar.pdb'), 'vhh')
cmd.load(str(STRUCT_DIR / 'vh_exemplar_aligned.pdb'), 'vh')

for obj, main_color, chains in [('vhh', 'vhh_blue', ['A']), ('vh', 'vh_purple', ['H'])]:
    cmd.hide('everything', obj)
    cmd.show('cartoon', f'{obj} and chain {"+".join(chains)}')
    cmd.color(main_color, f'{obj} and chain {"+".join(chains)}')
    sel_h = '+'.join(chains)
    cmd.color('iface_orange', f'{obj} and chain {sel_h} and resi {"+".join(map(str, SHIELDABLE_15))}')
    cmd.color('cdr3_red', f'{obj} and chain {sel_h} and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.show('spheres', f'{obj} and chain {sel_h} and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.color('tetrad_magenta', f'{obj} and chain {sel_h} and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.set('sphere_scale', 0.6, obj)

cmd.show('cartoon', 'vh and chain L')
cmd.color('vl_pink', 'vh and chain L')
cmd.set('cartoon_transparency', 0.7, 'vh and chain L')

cmd.orient('vhh')
cmd.zoom('vhh', buffer=3)
view = cmd.get_view()

cmd.disable('vh')
cmd.enable('vhh')
cmd.set_view(view)
ray_png(OUT_DIR / 'Figure4A_panelA_vhh.png')

cmd.disable('vhh')
cmd.enable('vh')
cmd.set_view(view)
ray_png(OUT_DIR / 'Figure4A_panelA_vh.png')

print('Figure 4A Panel A rendered')
