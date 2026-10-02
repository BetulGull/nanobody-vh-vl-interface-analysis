"""
Figure 5 revision -- render Panels F/G (median VHH/VH surfaces) and H (median
largest-patch overlay), using EXACTLY the same rendering logic, color ramp,
and camera-sharing convention as the existing Panels C/D/E
(scripts/stage12_render_raybould_style_physchem_surface.py), applied to the
deterministic nearest-to-median exemplars (Figure5_exemplars.csv).
"""
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/figure5_median_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

cmd.set_color("rbblue", [0.02, 0.02, 0.55])
cmd.set_color("rbcyan", [0.0, 0.55, 0.75])
cmd.set_color("rbteal", [0.0, 0.65, 0.55])
cmd.set_color("rbgreen", [0.30, 0.85, 0.15])

cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.set("orthoscopic", 1)
cmd.set("ray_trace_mode", 0)
cmd.set("surface_quality", 1)
cmd.set("antialias", 2)
cmd.set("specular", 0.15)
cmd.set("ambient", 0.4)
cmd.set("ray_shadows", 0)


def color_by_hydropathy(obj_sel):
    cmd.color("white", obj_sel)
    exposed_interface_sel = f"({obj_sel}) and b>0.99 and b<2.01"
    cmd.spectrum("b", "rbblue_rbcyan_rbteal_rbgreen", exposed_interface_sel, minimum=1.0, maximum=2.0)


# Panel F: VHH median surface
cmd.delete("all")
cmd.load(str(STRUCT_DIR / "vhh_median.pdb"), "mol")
cmd.show("surface", "chain A")
color_by_hydropathy("chain A")
cmd.orient("chain A")
cmd.zoom("chain A", buffer=3)
view_f = cmd.get_view()
cmd.ray(1400, 1400)
cmd.png(str(OUT_DIR / "Figure5_panelF_vhh_median.png"), dpi=600)

# Panel G: VH median surface -- VL-free coordinates, same shared camera
cmd.delete("all")
cmd.load(str(STRUCT_DIR / "vh_median_only_aligned.pdb"), "mol")
cmd.show("surface", "chain H")
color_by_hydropathy("chain H")
cmd.set_view(view_f)
cmd.ray(1400, 1400)
cmd.png(str(OUT_DIR / "Figure5_panelG_vh_median.png"), dpi=600)

print("Panels F and G rendered with shared camera transform.")

# Panel H: largest hydrophobic patch overlay for the SAME median exemplars
# (residue sets validated exact-match against stage8_patch_topology.py + the
# frozen local_thr4.0_largest_patch_residue_count aggregate -- see
# figure5_exemplar_selection_audit.md)
largest_patch_vhh_median_resi = "110+111+112"
largest_patch_vh_median_resi = "52+71+115+118"

cmd.delete("all")
cmd.load(str(STRUCT_DIR / "vhh_median.pdb"), "mol")
cmd.show("surface", "chain A")
color_by_hydropathy("chain A")
cmd.color("black", f"chain A and resi {largest_patch_vhh_median_resi}")
cmd.set_view(view_f)
cmd.ray(1400, 1400)
cmd.png(str(OUT_DIR / "Figure5_panelH_vhh_median_patch.png"), dpi=600)

cmd.delete("all")
cmd.load(str(STRUCT_DIR / "vh_median_only_aligned.pdb"), "mol")
cmd.show("surface", "chain H")
color_by_hydropathy("chain H")
cmd.color("black", f"chain H and resi {largest_patch_vh_median_resi}")
cmd.set_view(view_f)
cmd.ray(1400, 1400)
cmd.png(str(OUT_DIR / "Figure5_panelH_vh_median_patch.png"), dpi=600)

print("Panel H (median largest-patch overlay) rendered.")
print("ALL MEDIAN PANELS RENDERED. Shared camera view:", view_f)
