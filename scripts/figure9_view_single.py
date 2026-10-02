"""View-only helper: interactively inspect ONE Figure9 3D module/group at a
time (env vhh-interface-viz). NOT a modification of
scripts/figure9_render_3d_modules.py -- duplicates its exact load/color/marker
logic (unchanged) for a single chosen module+group, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select with env var FIG9_WHICH (default: paneld_vhh):
  paneld_vhh / paneld_vh   -- structure-mapping panel (median exemplar, tetrad
                              gold spheres + hydropathy-shift orange / charge-shift
                              steel-blue position markers)
  panelfg_vhh / panelfg_vh -- PPC/PNC charge_class surface on the 28-position
                              former-interface set (blue=Arg/Lys/His, red=Asp/Glu,
                              white=other former-interface, gray=rest of domain)
  panelh_vhh / panelh_vh   -- APBS continuum electrostatic surface, +-5 kT/e
                              (red=negative, blue=positive)
"""
from pymol import cmd
from pathlib import Path
import csv
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
MED_DIR = ROOT / "data/interim/figure5_median_structures"
APBS_DIR = ROOT / "data/interim/stage8/apbs_cohort"
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure9"

VHH_BLUE = [0.19, 0.38, 0.80]
VH_PURPLE = [0.49, 0.31, 0.66]
TETRAD_GOLD = [0.54, 0.43, 0.0]
HYDRO_ORANGE = [0.91, 0.33, 0.18]
CHARGE_STEEL = [0.29, 0.44, 0.65]
POS_BLUE = [0.13, 0.33, 0.80]
NEG_RED = [0.80, 0.13, 0.13]
NEUTRAL_GRAY = [0.82, 0.82, 0.82]

for nm, rgb in [('vhh_blue', VHH_BLUE), ('vh_purple', VH_PURPLE), ('tetrad_gold', TETRAD_GOLD),
                 ('hydro_orange', HYDRO_ORANGE), ('charge_steel', CHARGE_STEEL),
                 ('pos_blue', POS_BLUE), ('neg_red', NEG_RED), ('neutral_gray', NEUTRAL_GRAY)]:
    cmd.set_color(nm, rgb)

TETRAD = ['42', '49', '50', '52']
FULL_28 = ['40', '42', '44', '49', '50', '52', '55', '66', '67', '68', '69', '103',
           '107', '108', '109', '110', '111', '111A', '111B', '111C', '112', '112A',
           '112B', '112C', '112D', '113', '118', '119', '120']


def sel(chain, positions):
    return f"chain {chain} and resi " + '+'.join(positions)


which = os.environ.get('FIG9_WHICH', 'paneld_vhh')
cmd.delete('all')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('orthoscopic', 1)
cmd.set('ray_shadows', 0)
cmd.set('ambient', 0.4)
cmd.set('specular', 0.15)
cmd.set('antialias', 2)
cmd.set('surface_quality', 1)

if which in ('paneld_vhh', 'paneld_vh'):
    with open(SRC / 'Figure9_position_tracks.csv') as f:
        track_rows = list(csv.DictReader(f))
    TRACKS = {}
    for r in track_rows:
        if r['in_top5'] == 'True':
            TRACKS.setdefault(r['track_name'], []).append(r['IMGT_position'])

    if which == 'paneld_vhh':
        obj, chain, color, fname = 'vhh_d', 'A', 'vhh_blue', 'vhh_median.pdb'
    else:
        obj, chain, color, fname = 'vh_d', 'H', 'vh_purple', 'vh_median_only_aligned.pdb'
    cmd.load(str(MED_DIR / fname), obj)
    cmd.hide('everything', obj)
    cmd.show('cartoon', obj)
    cmd.color(color, obj)
    cmd.set('cartoon_transparency', 0.35, obj)

    for pos in TETRAD:
        cmd.show('spheres', f'{obj} and chain {chain} and resi {pos} and name CA')
        cmd.set('sphere_scale', 0.9, f'{obj} and chain {chain} and resi {pos} and name CA')
        cmd.color('tetrad_gold', f'{obj} and chain {chain} and resi {pos} and name CA')

    hydro_pos = [p for p in TRACKS.get('hydropathy_shift', []) if p not in TETRAD]
    if hydro_pos:
        hsel = f'{obj} and chain {chain} and resi ' + '+'.join(hydro_pos) + ' and name CA'
        cmd.show('spheres', hsel)
        cmd.set('sphere_scale', 0.55, hsel)
        cmd.color('hydro_orange', hsel)

    charge_pos = [p for p in TRACKS.get('charge_shift', []) if p not in TETRAD]
    if charge_pos:
        csel = f'{obj} and chain {chain} and resi ' + '+'.join(charge_pos) + ' and name CA'
        cmd.show('spheres', csel)
        cmd.set('sphere_scale', 0.55, csel)
        cmd.color('charge_steel', csel)

    cmd.orient(obj)
    print(f'Figure9 Panel D ({obj}) -- gold=FR2 tetrad (42/49/50/52), '
          f'orange=hydropathy-shift hotspot, steel-blue=charge-shift hotspot.')

elif which in ('panelfg_vhh', 'panelfg_vh'):
    if which == 'panelfg_vhh':
        obj, chain, fname = 'vhh_c', 'A', 'vhh_median.pdb'
    else:
        obj, chain, fname = 'vh_c', 'H', 'vh_median_only_aligned.pdb'
    cmd.load(str(MED_DIR / fname), obj)
    POS_AA = ['ARG', 'LYS', 'HIS']
    NEG_AA = ['ASP', 'GLU']
    cmd.hide('everything', obj)
    cmd.show('surface', obj)
    cmd.color('neutral_gray', obj)
    present_28 = [p for p in FULL_28
                  if cmd.count_atoms(f'{obj} and chain {chain} and resi {p} and name CA') > 0]
    fi_sel = sel(chain, present_28)
    cmd.color('white', f'{obj} and {fi_sel}')
    for aa_list, color in [(POS_AA, 'pos_blue'), (NEG_AA, 'neg_red')]:
        aa_sel = ' or '.join([f'resn {aa}' for aa in aa_list])
        cmd.color(color, f'{obj} and {fi_sel} and ({aa_sel})')
    cmd.orient(obj)
    print(f'Figure9 Panel F/G ({obj}) -- charge_class surface: blue=Arg/Lys/His, '
          f'red=Asp/Glu, white=other former-interface position, gray=rest of domain.')

elif which in ('panelh_vhh', 'panelh_vh'):
    APBS_LIMIT = 5.0
    if which == 'panelh_vhh':
        label, instance = 'vhh_apbs', 'pdb_00006jb2-A'
    else:
        label, instance = 'vh_apbs', 'pdb_00003gi9-H-L'
    pdb_path = APBS_DIR / instance / f'{instance}.pdb'
    dx_path = APBS_DIR / instance / f'{instance}.pqr.dx'
    cmd.load(str(pdb_path), label)
    cmd.load(str(dx_path), f'{label}_map')
    cmd.hide('everything', label)
    cmd.show('surface', label)
    cmd.color('white', label)
    ramp_name = f'{label}_ramp'
    cmd.ramp_new(ramp_name, f'{label}_map', [-APBS_LIMIT, 0, APBS_LIMIT],
                 ['red', 'white', 'blue'])
    cmd.set('surface_color', ramp_name, label)
    cmd.orient(label)
    print(f'Figure9 Panel H ({label}) -- APBS electrostatic surface, '
          f'+-{APBS_LIMIT} kT/e (red=negative, blue=positive).')

else:
    raise SystemExit(f'Unknown FIG9_WHICH={which!r}')
