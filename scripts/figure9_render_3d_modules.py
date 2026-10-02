"""Figure 9 3D modules -- PyMOL script (env vhh-interface-viz).

Renders three groups of panels, all reusing existing frozen/derived data only:
  (D) Structure-mapping panel: median VHH/VH exemplars (reused from Figure 5) with
      hydropathy-shift / charge-shift / tetrad position markers (existing per-position
      values, ranked by a predeclared top-5 rule -- see Figure9_position_tracks.csv).
  (F/G 3D) PPC/PNC charge_class surface: per-residue charge_class is a fixed,
      amino-acid-identity lookup already defined in stage8_patch_scoring.py
      (ASP/GLU=-1.0, ARG/LYS=+1.0, HIS=+0.1) -- reused verbatim here for
      visualization only, restricted to the 28-position former-interface set.
      This is NOT a patch-membership decomposition (that data does not exist).
  (H 3D) APBS continuum surface: real APBS .dx potential grids (frozen Stage-8
      output) loaded and ramped onto the domain surface, shared symmetric
      color limits (+-5 kT/e, a standard APBS-PyMOL rendering convention,
      applied identically to both groups -- not a data-derived value).
"""
from pymol import cmd
from pathlib import Path
import csv

ROOT = Path("/Users/betulgul/Desktop/staj")
MED_DIR = ROOT / "data/interim/figure5_median_structures"
APBS_DIR = ROOT / "data/interim/stage8/apbs_cohort"
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure9"
OUT = ROOT / "results/final_manuscript_figures/main"

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


def base_settings():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('ambient', 0.4)
    cmd.set('specular', 0.15)
    cmd.set('antialias', 2)
    cmd.set('surface_quality', 1)


def ray_png(path, width=1400, height=1400, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


with open(SRC / 'Figure9_position_tracks.csv') as f:
    track_rows = list(csv.DictReader(f))

TRACKS = {}
for r in track_rows:
    if r['in_top5'] == 'True':
        TRACKS.setdefault(r['track_name'], []).append(r['IMGT_position'])

TETRAD = ['42', '49', '50', '52']

FULL_28 = ['40', '42', '44', '49', '50', '52', '55', '66', '67', '68', '69', '103',
           '107', '108', '109', '110', '111', '111A', '111B', '111C', '112', '112A',
           '112B', '112C', '112D', '113', '118', '119', '120']
# (28-entry set incl. CDR3 insertion codes as sequence-resolved in Stage-3/Figure-1
#  audit; loop residues beyond 112D up to 118 are also part of the 28-entry table,
#  restricted below to only positions actually present with CA atoms in each
#  exemplar's resolved coordinates.)


def sel(chain, positions):
    return f"chain {chain} and resi " + '+'.join(positions)


# ============================================================
# (D) Structure-mapping panel
# ============================================================
cmd.delete('all')
base_settings()
cmd.load(str(MED_DIR / 'vhh_median.pdb'), 'vhh_d')
cmd.load(str(MED_DIR / 'vh_median_only_aligned.pdb'), 'vh_d')

for obj, chain, color in [('vhh_d', 'A', 'vhh_blue'), ('vh_d', 'H', 'vh_purple')]:
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

cmd.orient('vhh_d')
view = cmd.get_view()
cmd.disable('vh_d'); cmd.enable('vhh_d'); cmd.set_view(view)
ray_png(OUT / 'Figure9_panelD_vhh_median_markers.png')
cmd.disable('vhh_d'); cmd.enable('vh_d'); cmd.set_view(view)
ray_png(OUT / 'Figure9_panelD_vh_median_markers.png')
print('Panel D structure mapping rendered')

# ============================================================
# (F/G) PPC/PNC charge_class surface on former-interface positions
# ============================================================
cmd.delete('all')
base_settings()
cmd.load(str(MED_DIR / 'vhh_median.pdb'), 'vhh_c')
cmd.load(str(MED_DIR / 'vh_median_only_aligned.pdb'), 'vh_c')

POS_AA = ['ARG', 'LYS', 'HIS']
NEG_AA = ['ASP', 'GLU']

for obj, chain, base_color in [('vhh_c', 'A', 'vhh_blue'), ('vh_c', 'H', 'vh_purple')]:
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

cmd.orient('vhh_c')
view = cmd.get_view()
cmd.disable('vh_c'); cmd.enable('vhh_c'); cmd.set_view(view)
ray_png(OUT / 'Figure9_panelFG_vhh_median_chargeclass.png')
cmd.disable('vhh_c'); cmd.enable('vh_c'); cmd.set_view(view)
ray_png(OUT / 'Figure9_panelFG_vh_median_chargeclass.png')
print('PPC/PNC charge_class surfaces rendered')

# ============================================================
# (H) APBS continuum electrostatic surface
# ============================================================
cmd.delete('all')
base_settings()

APBS_LIMIT = 5.0  # kT/e, standard APBS-PyMOL rendering convention (not data-derived)

vhh_apbs_id = 'pdb_00006jb2-A'
vh_apbs_id = 'pdb_00003gi9-H-L'

for label, instance, chain_sel in [('vhh_apbs', vhh_apbs_id, 'chain A'),
                                     ('vh_apbs', vh_apbs_id, 'all')]:
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
    ray_png(OUT / f'Figure9_panelH_{label}.png')

print('APBS continuum surfaces rendered (shared +-5 kT/e color limits)')
