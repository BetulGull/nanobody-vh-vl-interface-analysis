# 3D structure render targets (PyMOL, env: vhh-interface-viz)
#
# Usage:
#   make figure7     -> Figure7 panels B/C (VHH/VH complex, cartoon + antigen surface)
#   make figure9     -> Figure9 3D modules (panel D markers, F/G charge-class, H APBS)
#   make figure11    -> Figure11 6 exemplar structures (hotspot + tetrad markers)
#   make all         -> all three headless renders, in order
#   make check-env   -> just verify micromamba + the vhh-interface-viz env are present
#
# Interactive, WHOLE-SCRIPT (opens a real PyMOL window you can rotate/zoom by
# hand -- same scene, same colors, same camera setup as the headless targets,
# just not headless; every structure the script loads stays in the object
# panel, only toggled show/hide, so you can inspect any of them by clicking
# their name in the panel on the right):
#   make figure7-view      -> loads BOTH Figure7 complexes (vhh, vh)
#   make figure9-view      -> runs all 3 Figure9 modules; only the LAST one
#                             (Panel H APBS, both vhh_apbs+vh_apbs) is still
#                             loaded at the end (Panel D/F-G get replaced,
#                             matching the original script's own cmd.delete
#                             calls) -- use the figure9-view-* targets below
#                             to open Panel D or F/G interactively instead.
#   make figure11-view     -> loads all 6 Figure11 exemplars at once
#
# Interactive, SINGLE-STRUCTURE (opens ONLY the one structure/module named,
# pre-oriented, nothing else to toggle -- uses small view-only helper scripts
# in scripts/*_view_single.py that duplicate the exact color/selection logic
# of the frozen render scripts verbatim, just for one structure at a time and
# without the final cmd.ray()/cmd.png() save step):
#   make figure7-view-vhh            -> Figure7 panel B (VHH+antigen)
#   make figure7-view-vh             -> Figure7 panel C (VH+VL+antigen)
#   make figure9-view-paneld-vhh     -> Figure9 panel D, VHH median exemplar
#   make figure9-view-paneld-vh      -> Figure9 panel D, VH median exemplar
#   make figure9-view-panelfg-vhh    -> Figure9 panel F/G, VHH charge_class surface
#   make figure9-view-panelfg-vh     -> Figure9 panel F/G, VH charge_class surface
#   make figure9-view-panelh-vhh     -> Figure9 panel H, VHH APBS surface
#   make figure9-view-panelh-vh      -> Figure9 panel H, VH APBS surface
#   make figure11-view-vhh1          -> Figure11 pdb_00006hjy-G  (VHH)
#   make figure11-view-vhh2          -> Figure11 pdb_00006hjx-F  (VHH)
#   make figure11-view-vhh3          -> Figure11 pdb_00008dqu-A  (VHH)
#   make figure11-view-vh1           -> Figure11 pdb_00005ggv-H-L (VH)
#   make figure11-view-vh2           -> Figure11 pdb_00003wfd-H-L (VH)
#   make figure11-view-vh3           -> Figure11 pdb_00007l7d-H-L (VH)
#
# Figure1 / Figure4A / Figure4B / Figure5 / Figure6 (color schemes follow each
# figure's own caption):
#   make figure1-view-panela         -> Figure1 panel A state (VH purple, VL pink, teal CDR3)
#   make figure1-view-panelb         -> Figure1 panel B state (+ orange/yellow/magenta interface)
#   make figure4a-view-vhh           -> Figure4A VHH exemplar (red CDR3, orange 15-pos, magenta tetrad)
#   make figure4a-view-vh            -> Figure4A VH exemplar (same scheme, VL translucent)
#   make figure4b-view-vhh-short     -> Figure4B VHH, Q25 CDR3 length (gold/teal/green/red/magenta)
#   make figure4b-view-vhh-long      -> Figure4B VHH, Q75 CDR3 length
#   make figure4b-view-vh-short      -> Figure4B VH, Q25 CDR3 length (VL translucent)
#   make figure4b-view-vh-long       -> Figure4B VH, Q75 CDR3 length
#   make figure5-view-panelf-vhh     -> Figure5 panel F, VHH median hydropathy surface
#   make figure5-view-panelg-vh      -> Figure5 panel G, VH median hydropathy surface
#   make figure5-view-panelh-vhh     -> Figure5 panel H, VHH median + black largest-patch overlay
#   make figure5-view-panelh-vh      -> Figure5 panel H, VH median + black largest-patch overlay
#   make figure6-view-short          -> Figure6 VHH, Q25 CDR3 length (teal CDR3, orange 14-pos target)
#   make figure6-view-long           -> Figure6 VHH, Q75 CDR3 length
#
# Every interactive target BLOCKS the terminal until you close the PyMOL
# window (type "quit" in PyMOL's own command line at the bottom of its
# window, or just quit the PyMOL application / Cmd+Q). Rotate = left-click
# drag, zoom = scroll or right-click drag, pan = middle-click drag.
#
# Nothing here recomputes any scientific value -- these targets only re-run
# (or, for the single-structure targets, faithfully re-execute a trimmed copy
# of) the existing PyMOL rendering scripts in scripts/ against already-frozen
# inputs. The single-structure helper scripts are view-only and contain no new
# logic: every color, selection, and camera call in them is copied unchanged
# from the corresponding frozen script.

# Resolved once, at parse time: prefer micromamba on PATH if present, otherwise
# fall back to the known install location (this machine's ~/.local/bin/micromamba
# is NOT on the default Terminal PATH, unlike a login-shell dev environment).
MICROMAMBA := $(shell command -v micromamba 2>/dev/null)
ifeq ($(strip $(MICROMAMBA)),)
MICROMAMBA := $(HOME)/.local/bin/micromamba
endif

ENV        := vhh-interface-viz
PYMOL      := $(MICROMAMBA) run -n $(ENV) pymol -cq
PYMOL_VIEW := $(MICROMAMBA) run -n $(ENV) pymol -q

.PHONY: all figure7 figure9 figure11 figure7-view figure9-view figure11-view check-env \
        figure7-view-vhh figure7-view-vh \
        figure9-view-paneld-vhh figure9-view-paneld-vh \
        figure9-view-panelfg-vhh figure9-view-panelfg-vh \
        figure9-view-panelh-vhh figure9-view-panelh-vh \
        figure11-view-vhh1 figure11-view-vhh2 figure11-view-vhh3 \
        figure11-view-vh1 figure11-view-vh2 figure11-view-vh3 \
        figure1-view-panela figure1-view-panelb \
        figure4a-view-vhh figure4a-view-vh \
        figure4b-view-vhh-short figure4b-view-vhh-long figure4b-view-vh-short figure4b-view-vh-long \
        figure5-view-panelf-vhh figure5-view-panelg-vh figure5-view-panelh-vhh figure5-view-panelh-vh \
        figure6-view-short figure6-view-long

all: figure7 figure9 figure11

check-env:
	@test -x "$(MICROMAMBA)" || { \
		echo "ERROR: micromamba not found at $(MICROMAMBA) (and not on PATH)."; exit 1; }
	@"$(MICROMAMBA)" env list 2>/dev/null | grep -q "$(ENV)" || { \
		echo "ERROR: micromamba env '$(ENV)' not found (run: $(MICROMAMBA) env list)"; exit 1; }
	@echo "OK: micromamba ($(MICROMAMBA)) + $(ENV) env found."

figure7: check-env
	@test -f data/interim/figure7_structures/vhh_complex_antigen_aligned.pdb || { echo "ERROR: missing data/interim/figure7_structures/vhh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f data/interim/figure7_structures/vh_complex_antigen_aligned.pdb || { echo "ERROR: missing data/interim/figure7_structures/vh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_B_residue_classes.csv || { echo "ERROR: missing Figure7_structure_B_residue_classes.csv"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_C_residue_classes.csv || { echo "ERROR: missing Figure7_structure_C_residue_classes.csv"; exit 1; }
	$(PYMOL) scripts/figure7_render_revised.py
	@echo "Figure7 -> results/final_manuscript_figures/main/Figure7_panelB_vhh_complex.png, Figure7_panelC_vh_complex.png"

figure7-view: check-env
	@test -f data/interim/figure7_structures/vhh_complex_antigen_aligned.pdb || { echo "ERROR: missing data/interim/figure7_structures/vhh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f data/interim/figure7_structures/vh_complex_antigen_aligned.pdb || { echo "ERROR: missing data/interim/figure7_structures/vh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_B_residue_classes.csv || { echo "ERROR: missing Figure7_structure_B_residue_classes.csv"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_C_residue_classes.csv || { echo "ERROR: missing Figure7_structure_C_residue_classes.csv"; exit 1; }
	@echo "Opening interactive PyMOL window for Figure7 -- rotate/zoom with the mouse, quit the PyMOL window when done."
	$(PYMOL_VIEW) scripts/figure7_render_revised.py

figure9: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing data/interim/figure5_median_structures/vhh_median.pdb"; exit 1; }
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing data/interim/figure5_median_structures/vh_median_only_aligned.pdb"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pdb || { echo "ERROR: missing APBS VHH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pqr.dx || { echo "ERROR: missing APBS VHH exemplar .dx grid"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pdb || { echo "ERROR: missing APBS VH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pqr.dx || { echo "ERROR: missing APBS VH exemplar .dx grid"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure9/Figure9_position_tracks.csv || { echo "ERROR: missing Figure9_position_tracks.csv"; exit 1; }
	$(PYMOL) scripts/figure9_render_3d_modules.py
	@echo "Figure9 -> results/final_manuscript_figures/main/Figure9_panelD_*.png, Figure9_panelFG_*.png, Figure9_panelH_*.png"

figure9-view: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing data/interim/figure5_median_structures/vhh_median.pdb"; exit 1; }
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing data/interim/figure5_median_structures/vh_median_only_aligned.pdb"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pdb || { echo "ERROR: missing APBS VHH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pqr.dx || { echo "ERROR: missing APBS VHH exemplar .dx grid"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pdb || { echo "ERROR: missing APBS VH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pqr.dx || { echo "ERROR: missing APBS VH exemplar .dx grid"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure9/Figure9_position_tracks.csv || { echo "ERROR: missing Figure9_position_tracks.csv"; exit 1; }
	@echo "Opening interactive PyMOL window for Figure9 -- rotate/zoom with the mouse, quit the PyMOL window when done."
	$(PYMOL_VIEW) scripts/figure9_render_3d_modules.py

figure11: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	$(PYMOL) scripts/figure11_render_exemplars.py
	@echo "Figure11 -> results/final_manuscript_figures/main/Figure11_pdb_*.png (6 files)"

figure11-view: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening interactive PyMOL window for Figure11 -- rotate/zoom with the mouse, quit the PyMOL window when done."
	$(PYMOL_VIEW) scripts/figure11_render_exemplars.py

# ============================================================
# Single-structure interactive views (see header comment above)
# ============================================================

figure7-view-vhh: check-env
	@test -f data/interim/figure7_structures/vhh_complex_antigen_aligned.pdb || { echo "ERROR: missing vhh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_B_residue_classes.csv || { echo "ERROR: missing Figure7_structure_B_residue_classes.csv"; exit 1; }
	@echo "Opening Figure7 VHH complex (panel B) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG7_WHICH=vhh $(PYMOL_VIEW) scripts/figure7_view_single.py

figure7-view-vh: check-env
	@test -f data/interim/figure7_structures/vh_complex_antigen_aligned.pdb || { echo "ERROR: missing vh_complex_antigen_aligned.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure7/Figure7_structure_C_residue_classes.csv || { echo "ERROR: missing Figure7_structure_C_residue_classes.csv"; exit 1; }
	@echo "Opening Figure7 VH+VL complex (panel C) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG7_WHICH=vh $(PYMOL_VIEW) scripts/figure7_view_single.py

figure9-view-paneld-vhh: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing vhh_median.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure9/Figure9_position_tracks.csv || { echo "ERROR: missing Figure9_position_tracks.csv"; exit 1; }
	@echo "Opening Figure9 Panel D, VHH median exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=paneld_vhh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure9-view-paneld-vh: check-env
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing vh_median_only_aligned.pdb"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure9/Figure9_position_tracks.csv || { echo "ERROR: missing Figure9_position_tracks.csv"; exit 1; }
	@echo "Opening Figure9 Panel D, VH median exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=paneld_vh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure9-view-panelfg-vhh: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing vhh_median.pdb"; exit 1; }
	@echo "Opening Figure9 Panel F/G, VHH charge_class surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=panelfg_vhh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure9-view-panelfg-vh: check-env
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing vh_median_only_aligned.pdb"; exit 1; }
	@echo "Opening Figure9 Panel F/G, VH charge_class surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=panelfg_vh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure9-view-panelh-vhh: check-env
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pdb || { echo "ERROR: missing APBS VHH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00006jb2-A/pdb_00006jb2-A.pqr.dx || { echo "ERROR: missing APBS VHH exemplar .dx grid"; exit 1; }
	@echo "Opening Figure9 Panel H, VHH APBS surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=panelh_vhh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure9-view-panelh-vh: check-env
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pdb || { echo "ERROR: missing APBS VH exemplar PDB"; exit 1; }
	@test -f data/interim/stage8/apbs_cohort/pdb_00003gi9-H-L/pdb_00003gi9-H-L.pqr.dx || { echo "ERROR: missing APBS VH exemplar .dx grid"; exit 1; }
	@echo "Opening Figure9 Panel H, VH APBS surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG9_WHICH=panelh_vh $(PYMOL_VIEW) scripts/figure9_view_single.py

figure11-view-vhh1: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00006hjy-G (VHH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vhh1 $(PYMOL_VIEW) scripts/figure11_view_single.py

figure11-view-vhh2: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00006hjx-F (VHH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vhh2 $(PYMOL_VIEW) scripts/figure11_view_single.py

figure11-view-vhh3: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00008dqu-A (VHH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vhh3 $(PYMOL_VIEW) scripts/figure11_view_single.py

figure11-view-vh1: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00005ggv-H-L (VH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vh1 $(PYMOL_VIEW) scripts/figure11_view_single.py

figure11-view-vh2: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00003wfd-H-L (VH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vh2 $(PYMOL_VIEW) scripts/figure11_view_single.py

figure11-view-vh3: check-env
	@test -d data/interim/figure11_structures || { echo "ERROR: missing data/interim/figure11_structures/"; exit 1; }
	@test -f results/final_manuscript_figures/source_data/Figure11/Figure11_position_mapping.csv || { echo "ERROR: missing Figure11_position_mapping.csv"; exit 1; }
	@echo "Opening Figure11 exemplar pdb_00007l7d-H-L (VH) only -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG11_WHICH=vh3 $(PYMOL_VIEW) scripts/figure11_view_single.py

# ============================================================
# Figure1 / Figure4A / Figure4B / Figure5 / Figure6 single-structure views
# ============================================================

figure1-view-panela: check-env
	@test -f data/interim/final_manuscript_structures/figure1_native_fv.pdb || { echo "ERROR: missing figure1_native_fv.pdb"; exit 1; }
	@echo "Opening Figure1 panel A state -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG1_WHICH=panela $(PYMOL_VIEW) scripts/figure1_view_single.py

figure1-view-panelb: check-env
	@test -f data/interim/final_manuscript_structures/figure1_native_fv.pdb || { echo "ERROR: missing figure1_native_fv.pdb"; exit 1; }
	@echo "Opening Figure1 panel B state -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG1_WHICH=panelb $(PYMOL_VIEW) scripts/figure1_view_single.py

figure4a-view-vhh: check-env
	@test -f data/interim/figure4a_structures/vhh_exemplar.pdb || { echo "ERROR: missing vhh_exemplar.pdb"; exit 1; }
	@echo "Opening Figure4A VHH exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4A_WHICH=vhh $(PYMOL_VIEW) scripts/figure4a_view_single.py

figure4a-view-vh: check-env
	@test -f data/interim/figure4a_structures/vh_exemplar_aligned.pdb || { echo "ERROR: missing vh_exemplar_aligned.pdb"; exit 1; }
	@echo "Opening Figure4A VH exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4A_WHICH=vh $(PYMOL_VIEW) scripts/figure4a_view_single.py

figure4b-view-vhh-short: check-env
	@test -f data/interim/figure4b_structures/vhh_short.pdb || { echo "ERROR: missing vhh_short.pdb"; exit 1; }
	@echo "Opening Figure4B VHH short-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4B_WHICH=vhh_short $(PYMOL_VIEW) scripts/figure4b_view_single.py

figure4b-view-vhh-long: check-env
	@test -f data/interim/figure4b_structures/vhh_long.pdb || { echo "ERROR: missing vhh_long.pdb"; exit 1; }
	@echo "Opening Figure4B VHH long-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4B_WHICH=vhh_long $(PYMOL_VIEW) scripts/figure4b_view_single.py

figure4b-view-vh-short: check-env
	@test -f data/interim/figure4b_structures/vh_short.pdb || { echo "ERROR: missing vh_short.pdb"; exit 1; }
	@echo "Opening Figure4B VH short-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4B_WHICH=vh_short $(PYMOL_VIEW) scripts/figure4b_view_single.py

figure4b-view-vh-long: check-env
	@test -f data/interim/figure4b_structures/vh_long.pdb || { echo "ERROR: missing vh_long.pdb"; exit 1; }
	@echo "Opening Figure4B VH long-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG4B_WHICH=vh_long $(PYMOL_VIEW) scripts/figure4b_view_single.py

figure5-view-panelf-vhh: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing vhh_median.pdb"; exit 1; }
	@echo "Opening Figure5 Panel F, VHH median surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG5_WHICH=panelf_vhh $(PYMOL_VIEW) scripts/figure5_view_single.py

figure5-view-panelg-vh: check-env
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing vh_median_only_aligned.pdb"; exit 1; }
	@echo "Opening Figure5 Panel G, VH median surface -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG5_WHICH=panelg_vh $(PYMOL_VIEW) scripts/figure5_view_single.py

figure5-view-panelh-vhh: check-env
	@test -f data/interim/figure5_median_structures/vhh_median.pdb || { echo "ERROR: missing vhh_median.pdb"; exit 1; }
	@echo "Opening Figure5 Panel H, VHH median + largest-patch overlay -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG5_WHICH=panelh_vhh $(PYMOL_VIEW) scripts/figure5_view_single.py

figure5-view-panelh-vh: check-env
	@test -f data/interim/figure5_median_structures/vh_median_only_aligned.pdb || { echo "ERROR: missing vh_median_only_aligned.pdb"; exit 1; }
	@echo "Opening Figure5 Panel H, VH median + largest-patch overlay -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG5_WHICH=panelh_vh $(PYMOL_VIEW) scripts/figure5_view_single.py

figure6-view-short: check-env
	@test -f data/interim/final_manuscript_structures/fig6_vhh_short.pdb || { echo "ERROR: missing fig6_vhh_short.pdb"; exit 1; }
	@echo "Opening Figure6 VHH short-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG6_WHICH=short $(PYMOL_VIEW) scripts/figure6_view_single.py

figure6-view-long: check-env
	@test -f data/interim/final_manuscript_structures/fig6_vhh_long_aligned.pdb || { echo "ERROR: missing fig6_vhh_long_aligned.pdb"; exit 1; }
	@echo "Opening Figure6 VHH long-CDR3 exemplar -- rotate/zoom with the mouse, quit the PyMOL window when done."
	FIG6_WHICH=long $(PYMOL_VIEW) scripts/figure6_view_single.py
