"""Figure 11 -- 3D localization of former-interface hotspots across real,
deterministically-selected, framework-aligned VHH and conventional-VH
exemplar structures. PyMOL script (env vhh-interface-viz).

Reuses the exact common-framework Kabsch alignment already computed in
this session (data/interim/figure11_structures/*.pdb, RMSD 0.29-1.40A
over 90 common FR CA atoms to a single VHH reference frame) -- no new
alignment logic, no averaged/artificial structure.

Feature accent colors (matching Figure 10):
  Sequence divergence      -> teal
  Hydropathy shift         -> orange
  Charge shift             -> steel blue
  Native exposure          -> gray
  CDR3 shielding (VHH only)-> purple
  Antigen engagement       -> red

Multi-feature positions get a segmented ring (one thin arc per feature)
around a neutral white/light-gray sphere -- never assigned to only one
feature.
"""
from pymol import cmd
from pathlib import Path
import csv
import math

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure11_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"
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


def base_settings():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('ambient', 0.45)
    cmd.set('specular', 0.1)
    cmd.set('antialias', 2)


def ray_png(path, width=1300, height=1300, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


with open(SRC / 'Figure11_position_mapping.csv') as f:
    mapping_rows = list(csv.DictReader(f))

STRUCTS = [
    ('pdb_00006hjy-G', 'VHH', 'vhh1_pdb_00006hjy-G.pdb', 'A'),
    ('pdb_00006hjx-F', 'VHH', 'vhh2_pdb_00006hjx-F.pdb', 'A'),
    ('pdb_00008dqu-A', 'VHH', 'vhh3_pdb_00008dqu-A.pdb', 'A'),
    ('pdb_00005ggv-H-L', 'VH', 'vh1_pdb_00005ggv-H-L.pdb', 'H'),
    ('pdb_00003wfd-H-L', 'VH', 'vh2_pdb_00003wfd-H-L.pdb', 'H'),
    ('pdb_00007l7d-H-L', 'VH', 'vh3_pdb_00007l7d-H-L.pdb', 'H'),
]

cmd.delete('all')
base_settings()

view = None
for inst, grp, fname, chain in STRUCTS:
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

    # former-interface: restrained thin semi-transparent surface patch
    fi_sel = f'{obj} and chain {chain} and resi ' + '+'.join(sorted(FULL_25_STR))
    cmd.show('surface', fi_sel)
    cmd.color('former_iface', fi_sel)
    cmd.set('transparency', 0.68, fi_sel)

    # Important-position markers (Option B): every important position gets
    # the SAME neutral gold "hotspot" sphere -- multi-feature membership is
    # not encoded on the structure itself (it would require 4-5 colors on
    # one small residue and become unreadable); instead the full per-position
    # feature membership is given once in the shared key beneath the figure.
    # Canonical tetrad positions get an extra magenta halo (a larger,
    # semi-transparent sphere drawn first, with the smaller opaque hotspot
    # sphere on top) so the ring reads as a thin halo, not a blob.
    rows_here = [r for r in mapping_rows if r['structure_id'] == inst]
    for r in rows_here:
        pos = r['IMGT_position']
        ca_sel = f'{obj} and chain {chain} and resi {pos} and name CA'
        if cmd.count_atoms(ca_sel) == 0:
            continue
        if pos in TETRAD:
            # opaque magenta sphere slightly larger than the gold hotspot sphere,
            # drawn first and centered on the same CA point -- reads as a clean
            # thin magenta rim around the gold marker (no transparency needed).
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

    if view is None:
        cmd.orient(f'{obj} and chain {chain}')
        cmd.turn('y', 15)
        view = cmd.get_view()
    cmd.disable(obj)

for inst, grp, fname, chain in STRUCTS:
    obj = inst.replace('-', '_').replace('.', '_')
    for o2, _, _, _ in STRUCTS:
        cmd.disable(o2.replace('-', '_').replace('.', '_'))
    cmd.enable(obj)
    cmd.set_view(view)
    safe_name = inst.replace('/', '_')
    ray_png(OUT_DIR / f'Figure11_{safe_name}.png')

print('Figure 11 exemplar renders complete')
