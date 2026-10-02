# Figure 5 Revision — Exemplar Selection Audit

## 1-2. Median calculation and nearest-structure selection (verified)

Source: `data/processed/stage8_surface_patch_descriptors.parquet` (frozen, `status=='OK'`, N=533 VHH / 1429 VH — matches expected Ns exactly), metric `former_interface_vicinity_PSH` (unchanged, same column used by the existing high-95th-percentile exemplars).

| group | median PSH | selected exemplar | selected PSH | |distance| |
|---|---|---|---|---|
| VHH | 76.031 | pdb_00006jb2-A | 76.038 | 0.0069 |
| VH | 84.590 | pdb_00004k2u-I-M | 84.655 | 0.0646 |

Selection rule (identical to the existing 95th-percentile rule, verified by exact reproduction — see Section 3-4): eligible = `primary_interface_position_coverage>=0.90` and resolution present; nearest-to-target by absolute PSH distance; ties broken by coverage (desc) then resolution (asc) then INSTANCE (lexical). 394 VHH / 743 VH structures eligible (identical eligible pools to the 95th-percentile selection, since the eligibility filter itself doesn't depend on the target value).

## 3-4. 95th-percentile targets and exemplars — re-verified, unchanged

Re-derived independently in this session directly from the frozen Stage-8 table and the exact same selection rule: VHH target=96.696, selected pdb_00008f6v-A (PSH=96.516, distance=0.180); VH target=113.508, selected pdb_00006wo3-H-L (PSH=113.561, distance=0.052). **Exact match** to the frozen `stage12_raybould_style_high_psh_exemplars.csv` (same instances, same distances, same eligible-pool sizes: 394/743) — confirms the existing 95th-percentile exemplars were already selected correctly; no replacement was made.

## 5. C/D correspond exactly to E

Panel E reuses the existing, unmodified `stage12_panel_E_vhh_patch.png` / `stage12_panel_E_vh_patch.png` renders, built from the same `vhh_exemplar.pdb` (pdb_00008f6v-A) / `vh_exemplar_only_aligned.pdb` (pdb_00006wo3-H-L) files as Panels C/D. No new structure was selected for Panel E.

## 6. F/G correspond exactly to H

Panel H was rendered in this revision from `vhh_median.pdb` (pdb_00006jb2-A) and `vh_median_only_aligned.pdb` (pdb_00004k2u-I-M) — the identical files used for Panels F/G (`scripts/figure5_render_median_panels.py`, shared camera view captured once from Panel F and reused for G/H, matching the established Panel-C-camera-reuse convention).

## 7. No visual preference used

Both median exemplars were selected purely by the predeclared nearest-to-median rule before any rendering occurred; no alternative structure was inspected or substituted after seeing the render.

## 8-10. Frozen values / hydropathy scale / patch definition unchanged

- Stage-8 `former_interface_vicinity_PSH` values: read-only, unmodified.
- Hydropathy scale: identical 4-anchor `rbblue`(1.00)/`rbcyan`(1.33)/`rbteal`(1.66)/`rbgreen`(2.00) ramp and B-factor encoding convention (`stage8_patch_components.parquet` `kyte_doolittle_norm`, sentinel -1.0 outside the interface/non-exposed) reused verbatim from `scripts/stage12_render_raybould_style_physchem_surface.py`.
- Largest-connected-hydrophobic-patch definition: `src/qc/stage8_patch_topology.py::compute_topology` (4.0 A primary threshold, exposed+hydrophobic vicinity-member residues) reused unmodified. **Validated exact match** against the frozen `local_thr4.0_largest_patch_residue_count` aggregate for all four exemplars before use: VHH high-95=9, VH high-95=10 (both match the existing hardcoded Panel-E residue lists exactly, confirming the re-implementation is correct), VHH median=5, VH median=4 (both match the frozen per-instance aggregate in `final_figures_model_input.parquet` exactly).

## Additional note: VH median coverage

The VH median exemplar (pdb_00004k2u-I-M) has `primary_interface_position_coverage=0.929` (below the 1.0 seen for the other three exemplars, but above the 0.90 eligibility floor used by the existing, unmodified selection rule) — this is expected and not a selection error; it is simply the nearest-to-median structure among the same eligible pool the 95th-percentile selection already used.

## Verdict

All 10 required checks pass. Median exemplars are genuinely new, deterministic, correctly derived; existing 95th-percentile exemplars, PSH values, hydropathy scale, and patch definition are all confirmed unchanged.
