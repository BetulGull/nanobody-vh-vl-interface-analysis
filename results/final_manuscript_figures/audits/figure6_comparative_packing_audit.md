# Figure 6 — Comparative Packing Audit

## 1-4. Identical definitions (verified)

VHH and VH both use `src/qc/stage7_bsa.py::compute_pairwise_bsa` **called verbatim, unmodified** — the only difference between the two runs is the input atom list (`get_variable_domain_atoms(instance, 'H')` for the respective group's instance). This guarantees (1) identical BSA definition, (2) identical packing-target residue set (`NONADJACENT_PACKING_TARGET`, 14 positions — confirmed no substitution with the 15-position `FULL_SHIELDABLE_TARGET`, which the function does not even default to), (3) identical SASA parameters (Lee-Richards, 1.4 Å probe, FreeSASA 2.2.1, same `calc_sasa` call), (4) identical CDR3 definition (IMGT 105-117, all real residues, same `select_cdr3_atoms` function).

**Direct validation**: recomputing the frozen VHH `nonadj_BSA_pair_total` values via this exact procedure reproduced them with max deviation 2.3e-13 (floating-point noise) across all 533 VHH structures — confirming the function is being invoked identically for both groups, not re-implemented.

## 5. Missing-residue handling (verified identical)

Both groups: if either atom subset (CDR3 or target) is empty for a given structure, `status='TARGET_NOT_EVALUABLE'` and the structure is excluded — never zero-filled. Verified: 1962/1962 structures (533 VHH + 1429 VH) returned `status='OK'`; zero exclusions were needed for this cohort.

## 6-7. Native conventional-VH coordinates preserved; no relaxation/remodeling

`get_variable_domain_atoms(instance, 'H')` reads directly from the native paired VH-VL experimental mmCIF file and extracts only the H-chain atoms (IMGT 1-128, altloc-resolved, no hydrogens) — no coordinate transformation, minimization, or remodeling of any kind is applied. VL is never loaded into this atom list (not deleted from the source file, not geometrically removed from native coordinates — simply not included in this specific computation's input), matching the task's explicit "omitted from the calculation, not physically isolated" requirement.

## 8. β3 direction correctly labeled

`group_bin = 1` for VHH, `0` for VH (verified in `src/qc/` interaction-fitting code and in `Figure6_interaction_model.csv`: `slope_VH = β1`, `slope_VHH = β1 + β3`). β3 = 27.23 (raw BSA model) is therefore unambiguously "VHH slope minus VH slope" — confirmed by direct check `slope_VHH - slope_VH == beta3_interaction` (exact, within floating point).

## 9. Common-support analysis validity

Common-support range `[5, 26]` was derived directly from the VHH cohort's own observed CDR3-length range (`min`, `max` of `d[d.group=='VHH']['CDR3_length']`), not chosen arbitrarily. Restricting VH to this same range (N_VH drops from 1429 to 1425 — only 4 VH structures fall outside VHH's observed range) and refitting gives β3=26.55 (vs 27.23 full-range) — same sign, ~2.5% magnitude change, same significance level. Robust.

## 10. Normalized-sensitivity formula validity

- `BSA_per_CDR3_residue = BSA / CDR3_length` — arithmetically direct, valid for all evaluable structures (CDR3_length > 0 always, by construction of the CDR3 definition).
- `packing_fraction_if_valid = BSA_CDR3_side / SASA_A_alone` — a valid, physically meaningful dimensionless fraction: the numerator (CDR3-side burial) is a strict subset/component of the denominator (CDR3's own total isolated surface area), both computed by the same `compute_pairwise_bsa` call with no new SASA run required. Bounded in [0, 1] by construction (a residue subset's burial cannot exceed the whole surface it's drawn from) — verified numerically: `packing_fraction_if_valid` observed range fully within [0,1] for all 1962 structures.

## Nonlinearity note

A LOWESS-vs-linear-fit comparison showed a larger deviation for the VH group (max ~291 Å², vs ~105 Å² for VHH) — this is driven by the sparse, long tail of the VH CDR3-length distribution (up to 63 residues, vs VHH's max of 26), where very few structures exist to anchor a smooth fit. Rather than fitting a full spline sensitivity model, this concern is directly addressed by the common-support analysis (Section 9 above), which restricts to the well-populated, comparable range and shows the same conclusion — the primary linear interaction result is not an artifact of the sparse VH tail.

## Verdict

All 10 required checks pass. The conventional-VH BSA calculation is a faithful, exact-procedure extension of the frozen VHH definition, not a new or divergent metric.
