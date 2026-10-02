"""View-only helper: interactively inspect ONE Figure5 median-surface panel at
a time (env vhh-interface-viz). NOT a modification of
scripts/figure5_render_median_panels.py -- duplicates its exact
load/color/patch logic (unchanged) for a single chosen panel, skips the final
cmd.ray()/cmd.png() save step, and leaves the PyMOL GUI open so it can be
rotated/zoomed by hand.

Select with env var FIG5_WHICH (default: panelf_vhh):
  panelf_vhh -- Panel F: VHH median surface, hydropathy spectrum (blue-cyan-teal-green)
  panelg_vh  -- Panel G: VH median surface (VL-free), same spectrum
  panelh_vhh -- Panel H: VHH median, largest hydrophobic patch overlaid in black
  panelh_vh  -- Panel H: VH median, largest hydrophobic patch overlaid in black
"""
from pymol import cmd
from pathlib import Path
import os

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure5_median_structures"

cmd.set_color("rbblue", [0.02, 0.02, 0.55])
cmd.set_color("rbcyan", [0.0, 0.55, 0.75])
cmd.set_color("rbteal", [0.0, 0.65, 0.55])
cmd.set_color("rbgreen", [0.30, 0.85, 0.15])


def color_by_hydropathy(obj_sel):
    cmd.color("white", obj_sel)
    exposed_interface_sel = f"({obj_sel}) and b>0.99 and b<2.01"
    cmd.spectrum("b", "rbblue_rbcyan_rbteal_rbgreen", exposed_interface_sel, minimum=1.0, maximum=2.0)


LARGEST_PATCH_VHH = "110+111+112"
LARGEST_PATCH_VH = "52+71+115+118"

which = os.environ.get('FIG5_WHICH', 'panelf_vhh')

cmd.delete('all')
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.set("orthoscopic", 1)
cmd.set("ray_trace_mode", 0)
cmd.set("surface_quality", 1)
cmd.set("antialias", 2)
cmd.set("specular", 0.15)
cmd.set("ambient", 0.4)
cmd.set("ray_shadows", 0)

if which == 'panelf_vhh':
    cmd.load(str(STRUCT_DIR / "vhh_median.pdb"), "mol")
    cmd.show("surface", "chain A")
    color_by_hydropathy("chain A")
    cmd.orient("chain A")
    cmd.zoom("chain A", buffer=3)
    print('Figure5 Panel F -- VHH median surface, hydropathy spectrum (blue=hydrophilic .. green=hydrophobic).')
elif which == 'panelg_vh':
    cmd.load(str(STRUCT_DIR / "vh_median_only_aligned.pdb"), "mol")
    cmd.show("surface", "chain H")
    color_by_hydropathy("chain H")
    cmd.orient("chain H")
    cmd.zoom("chain H", buffer=3)
    print('Figure5 Panel G -- VH median surface (VL-free), hydropathy spectrum.')
elif which == 'panelh_vhh':
    cmd.load(str(STRUCT_DIR / "vhh_median.pdb"), "mol")
    cmd.show("surface", "chain A")
    color_by_hydropathy("chain A")
    cmd.color("black", f"chain A and resi {LARGEST_PATCH_VHH}")
    cmd.orient("chain A")
    cmd.zoom("chain A", buffer=3)
    print('Figure5 Panel H -- VHH median, largest hydrophobic patch (resi 110/111/112) in black.')
elif which == 'panelh_vh':
    cmd.load(str(STRUCT_DIR / "vh_median_only_aligned.pdb"), "mol")
    cmd.show("surface", "chain H")
    color_by_hydropathy("chain H")
    cmd.color("black", f"chain H and resi {LARGEST_PATCH_VH}")
    cmd.orient("chain H")
    cmd.zoom("chain H", buffer=3)
    print('Figure5 Panel H -- VH median, largest hydrophobic patch (resi 52/71/115/118) in black.')
else:
    raise SystemExit(f'Unknown FIG5_WHICH={which!r}')
