"""
Stage 12 — Raybould-style molecular surface rendering (Panels C/D) plus
structural reference (Panel B) and patch-outline overlay (Panel E).

Run in the vhh-interface-viz environment (PyMOL available there):
  micromamba run -n vhh-interface-viz python3 scripts/stage12_render_raybould_style_physchem_surface.py

Colour logic (adapted from Raybould et al. 2019 Fig 2B, PMC6410772):
  - Surface OUTSIDE the frozen 28-position empirical former VH-VL interface: white.
  - Surface INSIDE the interface AND surface-exposed (Stage-8 frozen classification):
    coloured by the frozen Stage-8 normalized Kyte-Doolittle hydropathy value
    (range [1,2]), using a 4-anchor blue->cyan->teal->green ramp reconstructed
    for presentation (exact RGB values not given in the source text -- see
    stage12_raybould_fig2_visualization_source_audit.md Section J-notes).
  - Interface residues NOT surface-exposed (buried): also left white (color
    represents EXPOSED hydropathy only, matching the source paper's own logic
    -- Raybould's PSH itself is computed only over exposed residues).
Colour encodes per-residue hydropathy; it does NOT encode the PSH score itself
(a pairwise spatial quantity) -- see caption / Section 15 of the task spec.
"""
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/stage12_raybould_structures"
FIG_DIR = ROOT / "results/figures/stage12"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 4-anchor Raybould-style palette (graphical reconstruction, documented as such)
cmd.set_color("rbblue", [0.02, 0.02, 0.55])      # 1.00 hydrophilic
cmd.set_color("rbcyan", [0.0, 0.55, 0.75])       # 1.33
cmd.set_color("rbteal", [0.0, 0.65, 0.55])       # 1.66
cmd.set_color("rbgreen", [0.30, 0.85, 0.15])     # 2.00 hydrophobic

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
    """White for sentinel(-1)/buried/outside-interface; ramp for real [1,2] values."""
    cmd.color("white", obj_sel)
    exposed_interface_sel = f"({obj_sel}) and b>0.99 and b<2.01"
    cmd.spectrum("b", "rbblue_rbcyan_rbteal_rbgreen", exposed_interface_sel, minimum=1.0, maximum=2.0)


def render_surface_panel(pdb_path, chain_show, chain_hide, out_png, view_name=None):
    cmd.delete("all")
    cmd.load(str(pdb_path), "mol")
    if chain_hide:
        cmd.hide("everything", f"chain {chain_hide}")
        cmd.disable  # no-op, keep chain in session (not deleted) but not shown
    cmd.show("surface", f"chain {chain_show}")
    color_by_hydropathy(f"chain {chain_show}")
    cmd.set("transparency", 0.0)
    cmd.orient(f"chain {chain_show}")
    cmd.zoom(f"chain {chain_show}", buffer=3)
    if view_name:
        cmd.set_view(view_name)
    cmd.ray(1400, 1400)
    cmd.png(str(out_png), dpi=600)
    return cmd.get_view()


def main():
    # Panel C: VHH exemplar surface
    view_c = render_surface_panel(STRUCT_DIR / "vhh_exemplar.pdb", "A", None,
                                    FIG_DIR / "stage12_panel_C_vhh_surface.png")

    # Panel D: VH exemplar surface -- rendered from VH-ONLY coordinates (no VL atoms
    # present at all), matching this project's established VH_WITHOUT_VL_GEOMETRIC
    # surface-calculation convention (Stage 5/8: VL atoms are never included in the
    # SASA/surface calculation itself, not merely hidden after the fact -- loading VL
    # atoms into the same PyMOL object and only hiding them produces a spurious cavity
    # where VL used to pack, which is not the intended geometric counterfactual).
    # The VH+VL combined file remains on disk (vh_exemplar_with_vl_aligned.pdb) as the
    # full scientific source-data record (validation item 10: VL hidden, not deleted
    # from scientific source data) and is used for Panel B's structural reference.
    cmd.delete("all")
    cmd.load(str(STRUCT_DIR / "vh_exemplar_only_aligned.pdb"), "mol")
    cmd.show("surface", "chain H")
    color_by_hydropathy("chain H")
    cmd.set_view(view_c)  # SAME camera as panel C -- no independent rotation
    cmd.ray(1400, 1400)
    cmd.png(str(FIG_DIR / "stage12_panel_D_vh_surface.png"), dpi=600)

    print("Panels C and D rendered with shared camera transform.")

    # Panel B: structural reference -- native VH+VL, interface highlighted, no hydropathy coloring
    cmd.delete("all")
    cmd.load(str(STRUCT_DIR / "vh_exemplar_with_vl_aligned.pdb"), "ref")
    cmd.hide("everything")
    cmd.show("cartoon", "ref")
    cmd.color("gray70", "ref and chain H")
    cmd.color("lightpink", "ref and chain L")
    cmd.set("cartoon_transparency", 0.35, "ref and chain L")
    interface_resi = "40+42+44+49+50+52+55+66+67+68+69+103+107+108+109+110+111+112+113+114+115+116+118+119+120"
    cmd.color("orange", f"ref and chain H and resi {interface_resi}")
    cmd.set_view(view_c)
    cmd.ray(1400, 1400)
    cmd.png(str(FIG_DIR / "stage12_panel_B_structural_reference.png"), dpi=600)

    # Panel B inset: VHH alone, interface + CDR3 highlighted
    cmd.delete("all")
    cmd.load(str(STRUCT_DIR / "vhh_exemplar.pdb"), "vhhref")
    cmd.hide("everything")
    cmd.show("cartoon", "vhhref")
    cmd.color("gray70", "vhhref")
    cmd.color("orange", f"vhhref and resi {interface_resi}")
    cmd.color("cyan", "vhhref and resi 105-117")
    cmd.set_view(view_c)
    cmd.ray(900, 900)
    cmd.png(str(FIG_DIR / "stage12_panel_B_inset_vhh.png"), dpi=600)

    print("Panel B (structural reference) and inset rendered.")

    # Panel E: largest hydrophobic patch outline (VHH and VH), same coloring + dark contour residues
    largest_patch_vhh_resi = "2+105+108+110+112+113+116+118"
    largest_patch_vh_resi = "59+62+110+111+112+114+115+117+118"

    cmd.delete("all")
    cmd.load(str(STRUCT_DIR / "vhh_exemplar.pdb"), "mol")
    cmd.show("surface", "chain A")
    color_by_hydropathy("chain A")
    cmd.color("black", f"chain A and resi {largest_patch_vhh_resi}")
    cmd.set_view(view_c)
    cmd.ray(1400, 1400)
    cmd.png(str(FIG_DIR / "stage12_panel_E_vhh_patch.png"), dpi=600)

    cmd.delete("all")
    cmd.load(str(STRUCT_DIR / "vh_exemplar_only_aligned.pdb"), "mol")
    cmd.show("surface", "chain H")
    color_by_hydropathy("chain H")
    cmd.color("black", f"chain H and resi {largest_patch_vh_resi}")
    cmd.set_view(view_c)
    cmd.ray(1400, 1400)
    cmd.png(str(FIG_DIR / "stage12_panel_E_vh_patch.png"), dpi=600)

    print("Panel E (largest patch outline) rendered.")
    print("ALL PANELS RENDERED with shared camera view:", view_c)


if __name__ == '__main__':
    main()
