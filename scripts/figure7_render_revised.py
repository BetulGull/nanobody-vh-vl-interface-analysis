"""Figure 7 revision -- explicit 28-position former-interface engagement coloring.
PyMOL script (env vhh-interface-viz). Antigen-based alignment reused from the
matched-epitope validation gate (RMSD=2.34A, N=176 common antigen CA).

RENDERING FIX (this pass): restore the whole-complex cartoon representation
from the original pre-revision Figure 7 script (scripts/final_manuscript_
render_structures.py::figure7) -- antibody domains as CARTOON (not opaque
surface), antigen as SURFACE, whole-complex camera (cmd.orient/zoom on the
entire complex object, buffer=5) -- while keeping the NEW color semantics
(gray antigen, dark-gray matched epitope, yellow/red/orange classified
antibody residues as LOCAL surface patches on top of the cartoon). No
scientific values, contacts, classifications, or alignment are recomputed.

Color logic (mutually exclusive, F=28-pos former-interface, A=antigen-contacting):
  F ∩ A -> red (former-interface, antigen-engaged)
  F \\ A -> yellow (former-interface, not engaged in this structure)
  A \\ F -> orange (other antibody-antigen contact)
  matched epitope (antigen)  -> dark gray
  antigen (rest)             -> light gray
"""
from pymol import cmd
from pathlib import Path
import csv

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure7_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"
SRC = ROOT / "results/final_manuscript_figures/source_data/Figure7"

COL = {
    'vhh_blue': [0.19, 0.38, 0.80], 'vh_purple': [0.49, 0.31, 0.66], 'vl_pink': [0.95, 0.65, 0.77],
    'antigen_gray': [0.78, 0.78, 0.78], 'epitope_darkgray': [0.44, 0.44, 0.44],
    'fi_yellow': [1.0, 0.83, 0.0], 'fi_red': [0.82, 0.08, 0.08], 'other_orange': [1.0, 0.60, 0.0],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)

MATCHED_EPITOPE_RESI = [41, 42, 43, 44, 45, 47, 66, 67, 69, 71, 72]  # antigen numbering, shared frame


def base_settings():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('ambient', 0.4)
    cmd.set('specular', 0.15)
    cmd.set('antialias', 2)
    cmd.set('cartoon_transparency', 0.0)
    cmd.set('surface_quality', 1)


def ray_png(path, width=1500, height=1500, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


def apply_residue_classes(obj, chain, csv_path):
    """Show ONLY the classified residues as local surface patches on top of
    the cartoon -- never the whole chain."""
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    for cls, color in [('former_interface_nonengaged', 'fi_yellow'),
                        ('former_interface_engaged', 'fi_red'),
                        ('other_antigen_contact', 'other_orange')]:
        positions = [r['IMGT_position'] for r in rows if r['render_class'] == cls]
        if positions:
            sel = '+'.join(str(p) for p in positions)
            patch_sel = f'{obj} and chain {chain} and resi {sel}'
            cmd.show('surface', patch_sel)
            cmd.color(color, patch_sel)


cmd.delete('all')
base_settings()
cmd.load(str(STRUCT_DIR / 'vhh_complex_antigen_aligned.pdb'), 'vhh')
cmd.load(str(STRUCT_DIR / 'vh_complex_antigen_aligned.pdb'), 'vh')

# -- VHH complex: antibody as CARTOON, antigen as SURFACE, classified
# residues as small surface patches overlaid on the cartoon.
cmd.hide('everything', 'vhh')
cmd.show('cartoon', 'vhh and chain A')
cmd.color('vhh_blue', 'vhh and chain A')
cmd.show('surface', 'vhh and chain G')
cmd.color('antigen_gray', 'vhh and chain G')
epi_sel = '+'.join(str(r) for r in MATCHED_EPITOPE_RESI)
cmd.color('epitope_darkgray', f'vhh and chain G and resi {epi_sel}')
apply_residue_classes('vhh', 'A', SRC / 'Figure7_structure_B_residue_classes.csv')

# -- VH complex: VH heavy chain + VL both as CARTOON (distinct colors),
# antigen as SURFACE, classified VH residues as small surface patches.
cmd.hide('everything', 'vh')
cmd.show('cartoon', 'vh and chain H')
cmd.color('vh_purple', 'vh and chain H')
cmd.show('cartoon', 'vh and chain L')
cmd.color('vl_pink', 'vh and chain L')
cmd.show('surface', 'vh and chain G')
cmd.color('antigen_gray', 'vh and chain G')
cmd.color('epitope_darkgray', f'vh and chain G and resi {epi_sel}')
apply_residue_classes('vh', 'H', SRC / 'Figure7_structure_C_residue_classes.csv')

# -- Whole-complex camera, restored from the original pre-revision Figure 7
# script (final_manuscript_render_structures.py::figure7): orient/zoom on
# the ENTIRE complex object (full antibody + full antigen), not a tight
# interface-only crop. Camera computed once on the VHH complex and reused
# for the VH complex, since both structures are already antigen-aligned
# into a common frame (Kabsch, RMSD=2.34A, N=176 -- reused, not recomputed).
cmd.orient('vhh')
cmd.zoom('vhh', buffer=5)
view = cmd.get_view()

cmd.disable('vh')
cmd.enable('vhh')
cmd.set_view(view)
ray_png(OUT_DIR / 'Figure7_panelB_vhh_complex.png')

cmd.disable('vhh')
cmd.enable('vh')
cmd.set_view(view)
ray_png(OUT_DIR / 'Figure7_panelC_vh_complex.png')

print('Figure 7 panels B/C rendered: cartoon antibody domains, surface antigen, '
      'localized yellow/red/orange residue-class surface patches, whole-complex camera')
