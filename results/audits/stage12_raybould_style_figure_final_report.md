# Final Paper Visualization Task — Raybould-Style 3D Physicochemical Surface Figure

## A. Exact Raybould visual method reproduced/adapted

Reproduced exactly (independently confirmed from the primary source text, PMC6410772, not from memory): the ≥7.5% side-chain relative-ASA exposure threshold (Shrake-Rupley), the Kyte-Doolittle hydrophobicity scale normalized to [1,2], the <7.5 Å PSH pairwise cutoff, the 4 Å vicinity radius, and the core visual grammar — physicochemical coloring restricted to the target region, **explicit white everywhere else** ("Heavy- and light-chain surfaces outside the CDR vicinity are colored in white" — direct quote). See `results/audits/stage12_raybould_fig2_visualization_source_audit.md` for the full source audit.

## B. Exact difference between Raybould target region and ours

Raybould: **CDR vicinity** (every surface-exposed IMGT CDR + anchor residue, plus other surface-exposed residues within 4 Å). This project: the frozen **Stage-3 28-position FULL_EMPIRICAL_INTERFACE** — never substituted with the tetrad, the 15-position shieldable subset, or the Stage-7 14-position packing target.

## C. Full cohort Ns

Panel A uses the complete frozen Stage-8 cohort: **VHH N=533, conventional VH N=1429** (all `status=='OK'`, verified against `stage8_surface_patch_descriptors.parquet` — no subsampling).

## D. PSH distributions used

`former_interface_vicinity_PSH` (Stage 8, frozen, unmodified) for both groups, plotted as proportion-normalized overlapping histograms with identical fixed bins (31 bins spanning the full combined range), matching the conceptual presentation of Raybould Fig. 2A.

## E. VHH 95th-percentile PSH and selected structure

Group 95th percentile = **96.696**. Selected structure: **pdb_00008f6v / pdb_00008f6v-A** (PSH=96.516, distance=0.180, coverage=1.0, resolution=2.30 Å). 394/533 VHH structures were eligible (coverage≥0.90); selection made by minimum absolute distance to the frozen percentile, no ties encountered.

## F. VH 95th-percentile PSH and selected structure

Group 95th percentile = **113.508**. Selected structure: **pdb_00006wo3 / pdb_00006wo3-H-L** (PSH=113.561, distance=0.052, coverage=1.0, resolution=2.38 Å). 743/1429 VH structures were eligible; no ties encountered.

## G. Alignment positions and RMSD

87 common IMGT framework (FR1-FR4) Cα positions aligned via Kabsch superposition. **RMSD = 1.358 Å.** Positions used, VHH chain, VH chain, and VL chain are recorded in `results/structures/stage12_raybould_style_alignment_audit.csv`. CDR3, antigen, VL, and the empirical interface itself were explicitly excluded from the alignment target set.

## H. Exact hydropathy normalization

Kyte-Doolittle raw values, min-max normalized to **[1.00, 2.00]** across the 20 canonical amino acids — identical formula/values already frozen and validated in Stage 8 (`KD_NORM` in `stage8_patch_scoring.py`), independently re-confirmed against the Raybould primary source text in this task's Section A.

## I. Exact exposure rule

Side-chain relative ASA **≥7.5%** (Shrake-Rupley-family algorithm via the real `psa` binary, Stage 8's validated pipeline) — reused unchanged, never recalculated for this figure. Confirmed to be byte-for-byte the same threshold stated in the Raybould primary source text.

## J. Exact surface-rendering software/settings

**PyMOL** (available in the local `vhh-interface-viz` micromamba environment; confirmed functional via a headless smoke-test render before committing to this path). White background, orthographic projection (`orthoscopic=1`), `ray_shadows=0`, fixed `ambient`/`specular` settings, identical `cmd.set_view()` camera matrix applied to every panel after being captured once from Panel C, 1400×1400 px ray-traced renders at 600 dpi PNG export. Raybould's own text does not specify rendering software — this is this project's own choice, documented as such.

## K. Residue mapping validation

`results/structures/stage12_raybould_style_surface_residue_mapping.csv` — 246 rows (one per frozen 28-position-interface residue per structure), every value (exposure, relative ASA, normalized hydropathy, Stage-3 interface frequency, tetrad/CDR3/shieldable-15 flags) read directly from already-frozen Stage-3/8 tables. Zero values recomputed for this figure except patch-COMPONENT LABELING (see L).

## L. Largest-patch validation

Per-residue largest-hydrophobic-connected-component membership was **not persisted per-residue** in the original Stage-8 run (only aggregate counts were saved). Regenerated using the **identical, unchanged** frozen Stage-8 topology function (`stage8_patch_topology.py`, same 4.0 Å primary threshold, same exposed-hydrophobic-residue definition) for these 2 exemplars only. **Validated exact match against the frozen aggregate statistic**: VHH regenerated component size = 9 = frozen `local_thr4.0_largest_patch_residue_count`; VH regenerated = 10 = frozen value. This is a label-regeneration under unchanged parameters, not a new scientific calculation — documented explicitly in the validation table (check 19).

## M. Final figure paths

- `results/figures/stage12/stage12_fig_raybould_style_vhh_vs_vh_psh_surface.png` (300 dpi)
- `results/figures/stage12/stage12_fig_raybould_style_vhh_vs_vh_psh_surface.pdf`
- `results/figures/stage12/stage12_fig_raybould_style_vhh_vs_vh_psh_surface.tiff` (300 dpi)
- Individual panel PNGs (A-E, 600 dpi molecular surfaces) also retained in `results/figures/stage12/`.
- Rendering scripts: `scripts/stage12_render_raybould_style_physchem_surface.py` (PyMOL panels), `scripts/stage12_raybould_panel_A_and_assembly.py` (histogram + assembly), `src/qc/stage12_raybould_alignment.py`, `src/qc/stage12_raybould_write_structures.py`.

## N. Confirmation exemplars were selected BEFORE visualization

Confirmed: `results/structures/stage12_raybould_style_high_psh_exemplars.csv` was written and frozen (with `selection_rank`, `distance_to_95th_percentile`, and full eligibility/tie-break documentation) as a fully separate, prior step before any PyMOL rendering script was written or executed. The selection used only frozen Stage-4/8 metadata (coverage, resolution, PSH) — no visual inspection occurred at any point before the exemplars were fixed.

## O. Confirmation no exemplar was replaced manually

Confirmed: the exemplars rendered in the final figure (`pdb_00008f6v-A`, `pdb_00006wo3-H-L`) are identical to the frozen selection in `stage12_raybould_style_high_psh_exemplars.csv`. No alternate structure was substituted at any point, including after a real rendering bug (VL-cavity artifact, see below) was found and fixed — the fix corrected the *rendering method* (excluding VL atoms from the surface calculation, matching the project's own VH_WITHOUT_VL_GEOMETRIC convention), not the *choice of structure*.

**One real technical issue found and fixed during rendering** (documented transparently, not hidden): the first Panel D render, built from a VH+VL combined PyMOL object with VL merely hidden, produced a spurious surface cavity where VL used to pack (PyMOL's surface algorithm accounts for all atoms in an object regardless of visibility). Fixed by rendering Panel D/E from VH-only coordinates (no VL atoms present at all), consistent with this project's established VH_WITHOUT_VL_GEOMETRIC surface-calculation convention (Stage 5/8) — VL is excluded from the *calculation*, not merely hidden after the fact. The full VH+VL coordinate file remains on disk (`vh_exemplar_with_vl_aligned.pdb`) as the complete scientific source-data record and is used for Panel B's structural-reference cartoon (where the cavity artifact does not occur).

## P. Confirmation no Stage 2-12 scientific result changed

Verified directly: Stage-4 VHH N=533, Stage-3 paired primary N=1429, Stage-8 patch descriptors N=1962 — all unchanged. Full repository test suite: **331/331 PASSED** (no new tests were required for this figure-only task; the existing Stage-12 test suite already guards the frozen upstream values this figure reads from). No Stage 2-12 output file was modified; only new Stage-12 figure/audit/structure files were added.

---

**Status: Raybould-style figure generation COMPLETE.** All 20 programmatic validation checks pass (`results/audits/stage12_raybould_style_surface_validation.csv`). No new scientific stage was launched.
