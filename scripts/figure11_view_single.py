"""View-only helper: interactively inspect ONE Figure11 exemplar structure at a
time (env vhh-interface-viz). NOT a modification of
scripts/figure11_render_exemplars.py -- duplicates its exact load/color/marker
logic (unchanged) for a single chosen structure, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select structure with env var FIG11_WHICH (default: vhh1):
  vhh1 = pdb_00006hjy-G   (VHH, dist_from_centroid=0.690)
  vhh2 = pdb_00006hjx-F   (VHH, dist_from_centroid=0.838)
  vhh3 = pdb_00008dqu-A   (VHH, dist_from_centroid=0.854)
  vh1  = pdb_00005ggv-H-L (VH,  dist_from_centroid=0.593)
  vh2  = pdb_00003wfd-H-L (VH,  dist_from_centroid=0.615)
  vh3  = pdb_00007l7d-H-L (VH,  dist_from_centroid=0.681)
"""
from pymol import cmd
from pathlib import Path
import csv
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure11_structures"
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure11"

VHH_BLUE = [0.19, 0.38, 0.80]
VH_PURPLE = [0.49, 0.31, 0.66]
VL_PINK = [0.95, 0.65, 0.77]
FORMER_IFACE = [1.0, 0.90, 0.55]
HOTSPOT_GOLD = [0.72, 0.53, 0.04]
TETRAD_HALO = [0.85, 0.10, 0.65]

for name, rgb in [('vhh_blue', VHH_BLUE), ('vh_purple', VH_PURPLE), ('vl_pink', VL_PINK),
                    ('former_iface', FORMER_IFACE), ('hotspot_gold', HOTSPOT_GOLD),
                    ('tetrad_halo', TETRAD_HALO)]:
    cmd.set_color(name, rgb)

FULL_25 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110,
           111, 112, 113, 114, 115, 116, 118, 119, 120}
FULL_25_STR = set(str(x) for x in FULL_25) | {'111A', '112A', '112B'}
TETRAD = {'42', '49', '50', '52'}

STRUCTS = {
    'vhh1': ('pdb_00006hjy-G', 'VHH', 'vhh1_pdb_00006hjy-G.pdb', 'A'),
    'vhh2': ('pdb_00006hjx-F', 'VHH', 'vhh2_pdb_00006hjx-F.pdb', 'A'),
    'vhh3': ('pdb_00008dqu-A', 'VHH', 'vhh3_pdb_00008dqu-A.pdb', 'A'),
    'vh1':  ('pdb_00005ggv-H-L', 'VH', 'vh1_pdb_00005ggv-H-L.pdb', 'H'),
    'vh2':  ('pdb_00003wfd-H-L', 'VH', 'vh2_pdb_00003wfd-H-L.pdb', 'H'),
    'vh3':  ('pdb_00007l7d-H-L', 'VH', 'vh3_pdb_00007l7d-H-L.pdb', 'H'),
}

which = os.environ.get('FIG11_WHICH', 'vhh1')
inst, grp, fname, chain = STRUCTS[which]

with open(SRC / 'Figure11_position_mapping.csv') as f:
    mapping_rows = list(csv.DictReader(f))

cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.45)
cmd.set('specular', 0.1)
cmd.set('antialias', 2)

obj = inst.replace('-', '_').replace('.', '_')
cmd.load(str(STRUCT_DIR / fname), obj)
cmd.hide('everything', obj)
cmd.show('cartoon', f'{obj} and chain {chain}')
main_color = 'vhh_blue' if grp == 'VHH' else 'vh_purple'
cmd.color(main_color, f'{obj} and chain {chain}')
cmd.set('cartoon_transparency', 0.15, f'{obj} and chain {chain}')

if grp == 'VH':
    cmd.show('cartoon', f'{obj} and chain L')
    cmd.color('vl_pink', f'{obj} and chain L')
    cmd.set('cartoon_transparency', 0.55, f'{obj} and chain L')

fi_sel = f'{obj} and chain {chain} and resi ' + '+'.join(sorted(FULL_25_STR))
cmd.show('surface', fi_sel)
cmd.color('former_iface', fi_sel)
cmd.set('transparency', 0.68, fi_sel)

rows_here = [r for r in mapping_rows if r['structure_id'] == inst]
for r in rows_here:
    pos = r['IMGT_position']
    ca_sel = f'{obj} and chain {chain} and resi {pos} and name CA'
    if cmd.count_atoms(ca_sel) == 0:
        continue
    if pos in TETRAD:
        halo_name = f'halo_{obj}_{pos}'
        cmd.pseudoatom(halo_name, pos=cmd.get_atom_coords(ca_sel))
        cmd.hide('everything', halo_name)
        cmd.show('spheres', halo_name)
        cmd.color('tetrad_halo', halo_name)
        cmd.set('sphere_scale', 0.56, halo_name)
    cmd.show('spheres', ca_sel)
    cmd.color('hotspot_gold', ca_sel)
    cmd.set('sphere_scale', 0.40, ca_sel)
cmd.hide('labels')

cmd.orient(f'{obj} and chain {chain}')
cmd.turn('y', 15)

print(f'Figure11 single-structure view: {inst} ({grp}) -- pale-yellow=former-interface '
      f'surface patch, gold spheres=important hotspots, magenta halo=FR2 tetrad.')
