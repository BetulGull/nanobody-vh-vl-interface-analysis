# Figure 6 — Existing Data Lineage Audit

## VHH points

- **Source file**: `data/processed/final_figures_model_input.parquet`, column `nonadj_BSA_pair_total`.
- **N**: 533 (100% non-null for the VHH group).
- **Metric definition**: pairwise mutual buried surface area between CDR3 (IMGT 105-117, all real residues) and the frozen 14-position `NONADJACENT_PACKING_TARGET` ({40,42,44,49,50,52,55,66,67,68,69,103,119,120}), computed by `src/qc/stage7_bsa.py::compute_pairwise_bsa` (frozen, unmodified) — see [`figure6_bsa_definition.md`](figure6_bsa_definition.md) for the exact formula.
- **Residue target**: `NONADJACENT_PACKING_TARGET`, 14 positions — confirmed this is the 14-position set, not the 15-position `SHIELDABLE_NON_CDR3_INTERFACE` used elsewhere in the project.
- **SASA implementation**: `src/qc/stage3_sasa.py::calc_sasa` (FreeSASA, Lee-Richards algorithm, 1.4 Å probe, project-wide validated engine, same one used in Stage 3/5/6).
- **Coordinate treatment**: native experimental VHH heavy-chain coordinates (IMGT 1-128, altloc-resolved, no hydrogens); no VL present in a VHH structure.
- **Status**: fully validated, frozen, unchanged in this revision.

## VH points (as they existed before this revision)

- **Source file**: same `final_figures_model_input.parquet`, same column `nonadj_BSA_pair_total`.
- **Exact N**: **0**. Verified directly: `d[d.group=='VH']['nonadj_BSA_pair_total'].notna().sum() == 0` — every one of the 1429 VH rows for this column is `NaN`.
- **Whether actually calculated**: **no**. `stage7_bsa.py` was only ever invoked on the VHH representative cohort (`src/qc/stage7_full_cohort.py` loads `stage4_vhh_sequence_representatives.parquet` exclusively — confirmed in the Figure-4 lineage investigation, same finding applies here since it is the same Stage-7 module).
- **Whether validated**: not applicable — no values exist to validate.
- **Whether the existing legend/labeling was correct**: **no**. The prior Figure 6 Panel A plotted `d[d.group=='VH']['CDR3_length']` against `d[d.group=='VH']['nonadj_BSA_pair_total']` — an all-`NaN` y-series. matplotlib silently drops non-finite points from a scatter call, so the VH legend entry existed but **zero VH points were ever actually rendered**. This is a real, confirmed lineage gap, not a fabrication of VH data — no incorrect VH values were shown, but the figure visually implied a completed VHH-vs-VH comparison that did not exist.

## Which case applies (Section 3 classification)

**Case B**: only validated VHH BSA existed. The conventional-VH counterpart is calculated in this revision as a targeted extension, reusing the identical, unmodified `stage7_bsa.compute_pairwise_bsa` function (see [`figure6_bsa_definition.md`](figure6_bsa_definition.md) and [`figure6_comparative_packing_audit.md`](figure6_comparative_packing_audit.md) for the validation that this reuse is exact, not a re-implementation).
