# Figure 9 — Sequence → Structure → Surface Audit

Verifies the 19 required checks for the integrated Figure 9 rebuild. Every check was verified by direct load/computation in this session, not by memory.

## 1. Positions unchanged (28 empirical former-interface entries)

`Figure9_position_level_bridge.csv` and `Figure9_sequence_logo_VHH/VH.csv` carry exactly the same 28 IMGT-position entries as the frozen `figure1_empirical_interface_position_audit.csv` (25 unique integers + 111A/112A/112B). No position added, dropped, or renumbered. Confirmed identical order to Figure 1/2.

## 2. 25-vs-28 handling stated correctly

All 28 sequence-resolved entries are shown throughout (sequence logos, tracks, physicochemical matrix). None are silently reduced to 25. `Figure9_position_level_bridge.csv::sequence_resolved` is `True` for all 28 rows, consistent with the earlier project-wide finding that all 28 Stage-3 entries are sequence-resolved (the 25-vs-28 distinction only matters for structural/BSA-style analyses elsewhere in the manuscript, not for this sequence-to-surface figure).

## 3. Identical logo scaling (Panels A/B)

`Figure9_sequence_logo_VHH.csv`/`_VH.csv` are byte-identical reuses of `Figure2_panelA/B_*_sequence_logo.csv` (frequencies sum to 1.0 per position, X dropped and renormalized, same amino-acid color palette `figure1_sequence_logo_color_mapping.csv`). Both logos share the same y-axis (0-1.0) and the same 28-position x-axis ordering — directly comparable by eye, as in Figures 1/2.

## 4. JSD unsigned

`JSD` in the physicochemical matrix (Panel C) and the `seq_divergence` track reuse `stage4_position_chemistry_comparison.csv::JSD` verbatim — Jensen-Shannon divergence is inherently non-negative and unsigned; rendered on a sequential (not diverging) color scale, correctly reflecting this.

## 5. Hydropathy/charge direction correct

Panel C explicitly labels both ends of the diverging scales: ΔHydropathy is captioned "VHH more hydrophilic ← 0 → VHH more hydrophobic" (i.e., positive delta = VHH higher Kyte-Doolittle = more hydrophobic); ΔCharge is captioned "VHH more negative ← 0 → VHH more positive" (positive delta = VHH's mean sequence-charge score is higher, i.e. more positively charged). Verified against the raw sign convention used to build `delta_hydropathy` (`mean_KD_VHH - mean_KD_VH`) and `delta_charge` (VHH mean charge score minus VH mean charge score). Example: position 50 delta_charge=+0.901 (VHH more positive — a known E→R-family shift, confirmed against the Figure 2 enrichment table), position 49 delta_charge=-0.750 (VHH more negative).

## 6. Exposure not confused with hydropathy

`VHH_native_exposure`/`VH_native_exposure` (mean SASA, Å²) are kept in entirely separate rows/tracks from `delta_hydropathy` (mean Kyte-Doolittle) throughout — never combined into a single "surface hydrophobicity" value at the position level. The whole-domain SASA×KD combination (Zhang-style) exists only in Figure 4B as a distinct, already-labeled quantity and is not reused here.

## 7. Shielding is geometric, not chemical

`VHH_CDR3_shielding` (from `stage6_position_shielding_summary.csv::mean_delta_sasa_shield`) is a ΔSASA geometric masking quantity (Å² of surface area occluded by CDR3 presence), not a chemical/hydropathy value. It is reported only for VHH (a VHH-specific structural concept) and only for the 15-position `SHIELDABLE_NON_CDR3_INTERFACE` subset — CDR3-range positions (105-117 by IMGT numbering, appearing as NaN, hatched in Panel C) have no shielding value by construction, since CDR3 cannot geometrically shield itself.

## 8. PPC/PNC prevalence separated from magnitude

Both PPC and PNC are reported as two independent statistical tests, never one boxplot over all values including structural zeros: (a) prevalence — Fisher's exact test on patch-positive vs not (`Figure9_PPC/PNC_prevalence.csv`); (b) conditional magnitude — Mann-Whitney U + rank-biserial effect size among patch-positive structures only (`Figure9_PPC/PNC_conditional_magnitude.csv`). Panel F/G renders both as separate sub-panels (bar chart for prevalence, boxplot restricted to nonzero values for magnitude).

## 9. Patch-membership tracks use only real validated membership data

`Figure9_PPC_position_membership.csv` and `Figure9_PNC_position_membership.csv` are explicit `NOT AVAILABLE FROM VALIDATED OUTPUTS` placeholders (all-NaN value columns, `status` and `reason` columns populated). Confirmed by full-file review of `stage8_patch_topology.py`: connected-component patch membership is implemented only for the hydrophobic (PSH) score (`HYDROPHOBIC_AA`-restricted union-find); no equivalent charge-connectivity function exists anywhere in the codebase for PPC/PNC (a pairwise-sum score, not a connected-component patch). No position-level PPC/PNC membership track appears in Panels A/B, C, or the annotation provenance file — only the two-part cohort-level statistics and a visualization of the existing per-residue `charge_class` amino-acid-identity lookup (not a patch-membership decomposition) are shown.

## 10. APBS components separate

`Figure9_APBS_summary.csv` reports `interface_potential_positive_mean` and `interface_potential_negative_mean` as two independent rows/tests, never netted into one score. Panel H's boxplot renders them as two visually separated groups ("Positive component" / "Negative component"), each with its own group comparison (pos: p=3.53e-06; neg: p=5.40e-38).

## 11. APBS not residue-additive

The APBS 3D surfaces (Panel H) render the raw `.pqr.dx` continuum potential grid directly ramped onto the domain surface (`cmd.ramp_new` + `surface_color`) — a spatial visualization of the existing frozen grid, not a decomposition into per-residue contributions. No "APBS residue score" or "residue contribution logo" was computed or shown anywhere in this figure, per the explicit prohibition. `Figure9_annotation_provenance.csv` row for "APBS 3D continuum surface" explicitly states this.

## 12. Same APBS color limits both groups

Both APBS surfaces (VHH `pdb_00006jb2-A`, VH `pdb_00003gi9-H-L`) use the identical fixed ramp range ±5 kT/e (`cmd.ramp_new(..., [-5, 0, 5], ['red','white','blue'])`), applied identically, no independent autoscaling. This is stated as a rendering convention (not a data-derived value) in `Figure9_exemplar_selection.csv` and the caption.

## 13. Exemplar selection deterministic

All 3D exemplar choices are recorded in `Figure9_exemplar_selection.csv` with `target_value`, `selected_value`, and `distance_from_target` columns:
- PSH/structure-mapping/PPC/PNC modules: reuse the Figure-5 median-PSH exemplars (VHH `pdb_00006jb2-A`, distance=0.0069; VH `pdb_00004k2u-I-M`, distance=0.0646) — selected before any rendering, per the instruction not to select new exemplars unnecessarily.
- APBS module: VHH reuses the same exemplar (verified to have a valid `.pqr.dx` grid on disk and non-null APBS means, distance=0). VH required a **new, necessary** deterministic substitution — `pdb_00004k2u-I-M` has NaN APBS values (confirmed no `.pqr.dx` file exists on disk for it) — so `pdb_00003gi9-H-L` was selected as nearest-to-joint-median (z-scored positive+negative distance=0.0228) among the APBS-valid AND grid-available VH pool. This deviation is fully documented, not silent.

## 14. 3D orientations comparable

Panel D (structure mapping) and the PSH/PPC-PNC modules use the same median-exemplar PDB files already Kabsch-aligned in the Figure-5 revision (VH aligned onto VHH frame, RMSD=1.270 Å, N=90) — framework-based alignment, not phenotype-based. The APBS module orients each domain independently via `cmd.orient` (a reasonable default view per structure); because APBS is a continuum property compared only via shared numeric color limits and the separate 2-component cohort statistic (not via visual side-by-side residue correspondence), exact frame-matching was not required for that specific module and no comparability claim beyond "same color scale" is made for it.

## 15. No outcome-prediction added

No solubility, aggregation, thermal/conformational stability, folding free energy, binding affinity, KD/kon/koff, or biological-activity calculation was performed or plotted anywhere in this figure. No predictor for these properties was installed or run. Confirmed by full review of every script and CSV produced in this revision (`figure9_render_3d_modules.py` and all Python blocks): only existing frozen Stage-3/4/5/6/8/9 quantities and their direct, transparent aggregations/statistics/visualizations were computed.

## 16. No frozen metric overwritten

All frozen source files (`stage4_position_chemistry_comparison.csv`, `stage5_interface_surface_chemistry_residues.parquet`, `stage6_position_shielding_summary.csv`, `stage8_surface_patch_descriptors.parquet`, `final_figures_model_input.parquet`, `stage9_framework_antigen_contact_frequency.csv`) were opened read-only in every script in this revision — none were written to. `Figure9_position_level_bridge.csv` values were spot-checked to exactly match their frozen sources (e.g., position 50 JSD=0.814 matches `stage4_position_chemistry_comparison.csv` exactly).

## 17. Central claim discipline maintained

The caption (`Figure9_caption.txt`) uses only "is consistent with," "spatially corresponds to," and "co-localizes with" language for the APBS layer, and explicitly states no causal or single-residue-additive claim is made from APBS. No claim of solubility, aggregation-resistance, stability, or binding affinity appears anywhere in the caption or this audit.

## 18. Provenance completeness

`Figure9_annotation_provenance.csv` (30 rows) covers every annotation in every panel/module (A/B logos, 6 tracks, Panel C's 4 rows, Panel D's 5 marker/structure elements, PSH/PPC/PNC/APBS 2D+3D elements including the explicit "not available" residue-membership rows), each with `figure_panel, annotation, IMGT_position_if_applicable, metric, source_file, source_column, existing_or_derived, derivation_formula_if_derived`.

## 19. Cohort sizes match frozen expectations

Verified directly: PSH/PPC/PNC cohort — VHH N=533, VH N=1429 (matches `stage8_surface_patch_descriptors.parquet` group counts exactly). APBS-valid cohort — VHH N=521, VH N=1386 (matches `final_figures_model_input.parquet` non-null `interface_potential_positive_mean`/`interface_potential_negative_mean` counts exactly, as stated in the task's own placeholder expectation).

## Summary of the one necessary, fully-documented deviation from "reuse the Figure-5 exemplars everywhere"

The Figure-5 VH median exemplar (`pdb_00004k2u-I-M`) lacks a valid APBS `.pqr.dx` grid (APBS was not successfully run to completion for this structure — confirmed by directory listing and log inspection). Reuse was attempted first; only the APBS module substitutes a new deterministic VH exemplar (`pdb_00003gi9-H-L`), selected by the same nearest-to-median rule used everywhere else in this project. All other modules (structure mapping, PSH, PPC/PNC charge_class surfaces) continue to use the original Figure-5 VH median exemplar unchanged.

## Verdict

All 19 required checks pass.
