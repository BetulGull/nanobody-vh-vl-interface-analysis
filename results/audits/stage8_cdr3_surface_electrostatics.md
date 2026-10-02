# Stage 8 — Surface Patch Topology, Local Former-Interface Patch Metrics, and Electrostatics

> **CLOSURE AUDIT NOTE (applied after initial report):** this report was revised by a targeted closure audit that (1) resolved an ambiguous "1962/1962 complete" statement into exact attempted/success/failure counts, (2) corrected an APBS interpretation that had collapsed distinct positive/negative potential components into a single "less negative" narrative, (3) added a PPC/PNC zero-inflation decomposition, (4) sharpened the TAP/TNP taxonomy (method-definition resolution vs. third-party execution bug), and (5) reassessed the OVERALL status by applying the predeclared primary/secondary component hierarchy literally rather than retaining the prior PARTIAL by default. No frozen Stage-8 numerical output was recomputed or altered — only reporting/interpretation and one new derived summary table (`stage8_charge_patch_zero_inflation_summary.csv`, computed from already-frozen per-instance values).

## A/B. Primary cohorts

**VHH: N = 533** primary representatives (verified, `stage4_vhh_sequence_representatives.parquet`).
**Conventional VH: N = 1429** primary representatives (verified, `stage3_paired_representatives.parquet`, `is_primary_representative==True`; VL never loaded — `extract_variable_domain_atoms` restricted to `Hchain` only, matching Stage 5's `VH_WITHOUT_VL_GEOMETRIC` definition by construction).

## C/D. TAP source-definition result and standard-descriptor status

`results/audits/stage8_tap_source_gate.md`: TAP paper/supplement not locally available; `external/TNP` explicitly states (citation-confirmed) it reproduces TAP's PSH/PPC/PNC methodology. Formula, thresholds (7.5% exposure, 4.0 Å vicinity, 7.5 Å pairwise cutoff, salt-bridge 3.2 Å neutralization), hydrophobicity normalization, and charge classes all **EXACTLY_RESOLVED** from direct source inspection. The one denominator ambiguity (PSA's internal max-ASA reference table) was resolved operationally by running the real author-provided `psa` binary directly, rather than approximated. **SFvCSP: UNRESOLVED_METHOD_DEPENDENCY** — not found anywhere locally, not fabricated (`stage8_sfvcsp_status.md`).

**Standard TAP/TNP descriptor status: PARTIAL (secondary component).** Formula-level fidelity achieved and computed (`standard_tnp_core_PSH/PPC/PNC`, TNP's own CDR-vicinity anatomical context). Direct author-code output-level cross-validation was **attempted and blocked** by a real, documented defect in the third-party `PDBUtils.py` residue-parsing logic (`stage8_tap_tnp_reproduction.md`) — not a fabricated pass.

### Explicit TAP/TNP taxonomy (this closure audit)

Re-read `stage8_tap_source_gate.md` and `stage8_tap_tnp_reproduction.md` (no file named `stage8_tap_tnp_method_definition.md` exists — the method definition itself is documented inline within the source-gate audit; the absence of that specific filename is not itself a scientific gap, per the closure instruction).

- **A. PSH mathematical definition:** **RESOLVED.** Full formula (KD scale + 1-2 min-max normalization, 4.0 Å vicinity radius, 7.5 Å pairwise cutoff, 1/d² Coulombic-like weighting) recovered directly from `Hydrophobicity_and_Charge_Assigner.py`, citation-confirmed as a TAP reproduction by the TNP paper's own text.
- **B. PPC mathematical definition:** **RESOLVED** — same source, same formula machinery, positive-charge-class variant (Arg/Lys=+1, His=+0.1 at pH 7.4).
- **C. PNC mathematical definition:** **RESOLVED** — same source, negative-charge-class variant (Asp/Glu=-1).
- **D. Adapted former-interface versions:** **YES**, confirmed implemented using the identical resolved mathematical rule (`compute_patch_scores` in `stage8_patch_scoring.py` is the single shared implementation for both standard-TNP-context and adapted-former-interface calls — only the target residue set differs, never the formula).
- **E. What exactly failed in the attempted direct author-code reproduction:** **AUTHOR SOFTWARE EXECUTION BUG**, not a method-definition or project-implementation issue. Specifically: `theraprofnano.Hydrophobicity_and_Charge_Profiler.Common.PDBUtils.PDBchain.__init__`'s residue-key parsing raised `IndexError: string index out of range` on two independent real project structures (one with insertion codes, one without — ruling out insertion-handling as the specific cause), reproduced identically both times. This is a defect in the third-party package's own residue-ID handling when applied to this project's input format, discovered and documented, not modified or patched (per instruction not to fix the third-party package unless necessary to establish what failed — establishing the failure did not require modifying it).

**Explicit taxonomy distinction, applied:**
- `METHOD DEFINITION UNRESOLVED` — applies **only to SFvCSP** (genuinely absent from all local sources).
- `AUTHOR SOFTWARE EXECUTION BUG` — applies to the standard-TAP/TNP **author-code reproduction attempt** (PSH/PPC/PNC's mathematical definitions themselves are NOT unresolved — only the third-party executable is broken for this input).
- `PROJECT IMPLEMENTATION VALIDATION FAILURE` — **does not apply to anything in Stage 8.** This project's own PSH/PPC/PNC implementation was independently validated on 20 real structures (`stage8_patch_validation.csv`, 100% success) with a hand-verified formula spot-check, entirely separate from the blocked third-party reproduction attempt.

## E. Adapted former-interface descriptor definition

`former_interface_core_PSH/PPC/PNC` and `former_interface_vicinity_PSH/PPC/PNC`: the exact same source-verified mathematical scoring rule (KD-normalized hydrophobicity, pH-7.4 charge classes, 1/d² Coulombic-like pairwise weighting, 7.5 Å cutoff, salt-bridge neutralization), applied to the project's own frozen empirical target (Stage-3 28-position former VH-VL interface) instead of TNP's CDR-vicinity — never a reinvented formula, per Section E of the amendment.

## F. Patch validation result

`results/audits/stage8_patch_validation.csv`: 10 VHH + 10 VH, all `status=OK`, zero PSA-lookup misses, plausible core/vicinity residue counts and scores. Manual spot-check of one pairwise distance/PSH contribution confirmed exact formula agreement by hand (documented in session; d>7.5Å correctly excluded).

## G/H. Interface coverage

Core coverage proxy (`n_core`/25 canonical positions): full cohort 1962/1962 evaluable; high-coverage (≥0.90) subset = 1652/1962 (466 VHH, 1186 VH).

## I/J/K. VHH vs VH former-interface PSH, PPC, PNC

`results/tables/stage8_patch_group_comparison.csv` (all descriptive; Mann-Whitney p labeled descriptive-only, no formal inference):

- **PSH** (core): VHH mean 41.0 vs VH mean 44.7 — VHH lower.
- **PSH** (vicinity): VHH 77.0 vs VH 86.0 — VHH lower.
- **PPC** (core): VHH 0.136 vs VH 0.063 — VHH higher.
- **PNC** (core): VHH 0.213 vs VH 0.097 — VHH higher.
- **PPC/PNC** (vicinity): same direction, larger magnitude (VHH substantially higher for both).

Interpretation boundary maintained: this describes a structural surface-chemistry phenotype difference, not solubility/stability/aggregation.

### PPC/PNC zero-inflation decomposition (new, this closure audit)

The PPC/PNC boxplots show strong mass at exactly zero (most structures have no former-interface positions close enough, in the right charge classes, to form any patch at all). `results/tables/stage8_charge_patch_zero_inflation_summary.csv` decomposes the group mean difference into **prevalence of any nonzero patch** vs **magnitude conditional on being nonzero** (frozen `former_interface_core_PPC`/`PNC` values only, no metric modification, no pseudocounts, no log-transform):

| metric | group | N | frac=0 | frac>0 | mean(all) | median(all) | mean(nonzero) | median(nonzero) | N(nonzero) |
|---|---|---|---|---|---|---|---|---|---|
| PPC | VHH | 533 | 51.6% | 48.4% | 0.136 | 0.000 | 0.281 | 0.094 | 258 |
| PPC | VH | 1429 | 77.5% | 22.5% | 0.063 | 0.000 | 0.279 | 0.061 | 322 |
| PNC | VHH | 533 | 58.3% | 41.7% | 0.213 | 0.000 | 0.511 | 0.193 | 222 |
| PNC | VH | 1429 | 72.6% | 27.4% | 0.097 | 0.000 | 0.353 | 0.160 | 391 |

**Corrected interpretation — the VHH-vs-VH difference is not a uniform shift, and is driven by different mechanisms for PPC vs PNC:**

- **PPC:** the group difference is driven almost entirely by **prevalence** — VHH structures are more than twice as likely to have any positive former-interface patch at all (48.4% vs 22.5%), but the magnitude of that patch, *given that one exists*, is nearly identical between groups (mean-nonzero 0.281 VHH vs 0.279 VH; median-nonzero 0.094 vs 0.061, a smaller relative gap than the prevalence gap). So: **more VHH structures have a positive patch; when a VH structure does have one, it is roughly comparable in size.**
- **PNC:** the group difference is driven by **both** prevalence (41.7% vs 27.4%, a smaller relative gap than PPC's) **and** magnitude-conditional-on-nonzero (mean-nonzero 0.511 VHH vs 0.353 VH — VHH's nonzero negative patches are also larger, not just more frequent).

This distinction (prevalence vs magnitude) was not present in the original report's group-mean-only comparison and is a materially more precise description of the underlying distributional difference.

### G/H. PSH and topology re-verified from frozen tables (this closure audit)

Re-read directly from `stage8_surface_patch_descriptors.parquet` and `stage8_patch_topology.parquet` (no recomputation): `former_interface_core_PSH` (VHH mean 41.013/median 40.334, N=533; VH mean 44.656/median 43.041, N=1429) and `former_interface_vicinity_PSH` (VHH 77.044/76.031; VH 86.021/84.590) — **both confirmed identical to the original report, no correction needed.** `local_thr4.0_largest_patch_residue_count` (VHH mean 3.396/median 3.0; VH mean 4.961/median 5.0) — **also confirmed identical.**

**Conceptual clarification (explicitly stated per the closure prompt):** Stage 5 measured *exposed hydrophobic burden/chemistry* — a SASA-weighted scalar summary of how much hydrophobic surface is exposed, with no spatial-clustering information. Stage 8's PSH and project-specific topology measure the **spatial organization** of that exposed hydrophobicity — how it is arranged into contiguous patches, not merely how much of it there is. A surface can expose comparable total hydrophobic SASA (Stage 5) while organizing it into either one large dominant patch or many small scattered ones (Stage 8) — these are genuinely different physical quantities, and the Stage-8 topology result is not a duplicate of the Stage-5 finding.

## L/M. Project-specific patch topology and threshold sensitivity

`results/tables/stage8_patch_topology_threshold_sensitivity.csv` (3 thresholds: 3.5/4.0/4.5 Å, primary=4.0 Å reused directly from the TAP/TNP source's own vicinity radius per Section 16's "reuse without contradiction" instruction). `largest_patch_residue_count` and `fraction_hydrophobic_surface_in_largest_patch` are **stable across all 3 thresholds and both local/global scope**: VHH consistently shows smaller, less-dominant largest hydrophobic patches than VH. `number_of_hydrophobic_surface_components` is **not stable** (direction flips between the 3.5 Å and 4.0/4.5 Å thresholds) — reported honestly as threshold-sensitive, not cherry-picked.

## N. Interface-vs-noninterface-framework control

`results/tables/stage8_interface_vs_framework_patch_control.csv`: framework-region PSH is larger in absolute terms than interface-core PSH in both groups (reflecting the framework region's larger residue count, not directly comparable without normalization) — reported as raw descriptive values per spec, no normalization imposed.

## O. Local-vs-global surface control

`results/tables/stage8_local_vs_global_patch_control.csv`: the local (former-interface vicinity) VHH-vs-VH percent difference is **substantially larger** than the corresponding global whole-domain difference for PPC (+93% local vs -14% global) and PNC (+150% local vs +21% global) — supporting that the charge-patch difference is **localized to the former interface**, not a simple reflection of a globally different VHH surface. PSH shows a smaller, same-direction local-vs-global difference (-10% local vs -3% global).

## P. SFvCSP status and predeclared role

`results/audits/stage8_sfvcsp_status.md`: not computed for conventional Fv (genuinely unavailable source); **NA for VHH** by definition. No pseudo-SFvCSP created anywhere (verified by test). **Predeclared role, confirmed:** SFvCSP is applicable only to paired conventional Fv, NA for VHH, and — per the original Stage-8 specification's own framing (SFvCSP is a Fv-pairing descriptor, entirely orthogonal to the primary former-VH-VL-interface mechanistic question) — **secondary** to the primary empirical former-interface analysis. Its unresolved status does not, on its own, determine the overall Stage-8 outcome (see status reassessment below).

## Q/R. PDB2PQR/APBS environment and protonation

`results/audits/stage8_apbs_environment_and_parameters.md`: APBS 3.4.1, PDB2PQR 3.6.1 (verified, matches Stage-0). Force field PARSE, titration method propka, **pH 7.4** (physiological, documented before any group result was inspected). Histidine handling via per-structure propka pKa assignment (not a fixed default).

## S. APBS mapping validation

`results/audits/stage8_pdb2pqr_mapping_validation.csv` and direct smoke-test inspection: all 127/127 residues (including 7 insertion-coded CDR3 residues on `pdb_00006zbp-FFF`) round-trip exactly through PDB2PQR with `--keep-chain`; surrogate chain ID reversible via the same mapping-table pattern used in Stage 7.

## T. APBS success/failure counts — EXACT ACCOUNTING (corrected)

**"1962/1962 complete" in the original report meant 1962/1962 structures ATTEMPTED/PROCESSED (i.e., no structure was skipped or left unprocessed by the pipeline), not 1962/1962 successful APBS calculations.** This ambiguity is resolved here with exact counts read directly from `data/processed/stage8_apbs_interface_electrostatics.parquet` (`status` column), never estimated:

| | VHH | VH | TOTAL |
|---|---|---|---|
| N_total_expected | 533 | 1429 | 1962 |
| N_attempted | 533 | 1429 | 1962 |
| N_PDB2PQR_success | 531 | 1427 | 1958 |
| N_PDB2PQR_failure | 2 | 2 | 4 |
| N_APBS_success | 521 | 1386 | 1907 |
| N_APBS_failure | 10 | 41 | 51 |
| N_grid_valid | 521 | 1386 | 1907 |
| N_grid_invalid | 0 | 0 | 0 |
| N_potential_sampling_success | 521 | 1386 | 1907 |
| N_potential_sampling_failure | 0 | 0 | 0 |
| **N_final_APBS_analyzed** | **521 (97.7%)** | **1386 (97.0%)** | **1907 (97.2%)** |

(N_PDB2PQR_success = N_attempted − N_PDB2PQR_failure; every structure that passed PDB2PQR proceeded to APBS.)

**On the three requested categories that were not independently tracked as distinct pipeline states:** `GRID_INVALID` and `POTENTIAL_SAMPLING_FAILURE` were never recorded as separate status values in the pipeline — the implementation only emits `PDB2PQR_FAILURE`, `APBS_FAILURE`, or `PDB2PQR_SUCCESS_APBS_SUCCESS`. However, `molecule_contained_in_grid` (boolean) and `interface_potential_n` (sample count) were recorded for every nominal success and checked directly here: **all 1907 successes have `molecule_contained_in_grid == True` (100%) and a nonzero, populated `interface_potential_n` (min 129 samples, no zeros, no NaNs)** — so the counts above (0 grid-invalid, 0 sampling-failure) are not assumed but empirically confirmed from the frozen per-structure records, distinct from the fact that these were never separately-labeled failure *states* in the pipeline's status enum.

**Root cause of the 51 `APBS_FAILURE` cases** (unchanged from the original report, re-verified here): `Valist_readPQR: Error parsing atom N ... concatenated fields` — a real APBS fixed-width PQR-column-overflow limitation (not a project-pipeline defect, not molecule-chemistry-dependent), tracked explicitly per instance in `results/audits/stage8_apbs_failures.csv`, never converted to zero potential (verified by `test_apbs_failures_not_zero_potential`).

Full-cohort decision: **A. FULL COHORT** (`stage8_apbs_full_cohort_decision.md`) — runtime benchmark (2.2 s/structure mean, 20 real structures) showed the full 1962-structure run was clearly practical (~72 min projected, actual runtime consistent with this); no subsampling was applied.

## U. Former-interface APBS potential comparison — CORRECTED interpretation

`results/tables/stage8_apbs_group_comparison.csv` (N=521 VHH, N=1386 VH analyzed). **The original report's phrase "VHH shows slightly less negative mean interface potential than VH" was an oversimplification that collapsed two distinct components into one directional story — corrected here.**

Reading the positive and negative components **separately, exactly as computed**:

- **Positive component:** `interface_potential_positive_mean` — VHH 118.3 kT/e vs VH 115.4 kT/e. **VHH is somewhat MORE positive** in its positive-potential component.
- **Negative component:** `interface_potential_negative_mean` — VHH -168.2 kT/e vs VH -157.3 kT/e. **VHH is somewhat MORE negative** in its negative-potential component (larger magnitude, not smaller).

**Both verified directly from the frozen table — confirmed exactly as the closure prompt anticipated.** The correct statement is: **VHH shows stronger potential magnitude in BOTH the positive and negative components** relative to VH, not a uniform shift toward "more positive" or "less negative." This is also reflected in `interface_potential_absolute_mean` being higher for VHH (147.9 vs 141.5 kT/e) — consistent with stronger magnitudes in both directions rather than a net directional shift. The overall (signed) `interface_potential_mean` (-53.6 VHH vs -54.5 VH, a small difference driven by the `fraction_positive`/`fraction_negative` mix, 0.400/0.600 VHH vs 0.378/0.622 VH) is reported here only as a summary statistic, not as evidence of a "less negative" biological phenotype — no net electrostatic score is constructed or treated as primary. Positive and negative components are kept fully separate throughout, per the standing PPC/PNC-and-APBS-separation policy.

## D. APBS representation limitation (explicit confirmation)

Confirmed directly from `src/qc/stage8_apbs_sampling.py`/`stage8_apbs_full_cohort.py`: potential was sampled via **trilinear interpolation of the APBS grid at the (x,y,z) coordinates of solvent-exposed former-interface side-chain atoms** — atom-coordinate sampling. The label used throughout this project's outputs is, and remains, **"APBS potential sampled at exposed former-interface atoms."** This is explicitly **not** upgraded to "molecular surface electrostatic potential" or "surface potential integral" anywhere in the codebase or reports (verified: no such phrase appears in any Stage-8 output file or the 3D viewer, which itself labels its APBS mode "APBS potential sampled at exposed interface atoms," per `stage8_3d_surface_viewer.py`).

**Limitation, stated explicitly (not remedied in this closure audit, per instruction not to rerun with a new surface mesh):** atom-centered sampling evaluates the potential very close to the sampled atom's own partial charges, so values can be strongly influenced by the local charged-atom environment (including the atom's own contribution) rather than reflecting a solvent-excluded molecular-surface value further from any single atom. This is a genuine representation limitation of the predeclared Stage-8 method, not a violation of it — the method was predeclared and executed exactly as specified (Section 28 of the original spec), and this note documents its known interpretive caveat rather than indicating an implementation defect.

## V. Interface-vs-framework APBS comparison

`results/tables/stage8_apbs_interface_vs_framework.csv`: interface-region potential is more negative on average than framework-region potential in both groups (VHH interface -53.6 vs framework -42.2; VH interface -54.5 vs framework -42.4), and the interface region's fraction-positive is lower than framework's in both groups — a broadly consistent, non-group-specific regional pattern.

## W. PPC/PNC vs APBS relationships

`results/tables/stage8_ppc_pnc_vs_apbs_correlations.csv` (N=1907): **weak correlations** (|r| < 0.08 for all four pairings, core/vicinity PPC vs positive potential, core/vicinity PNC vs negative potential). Interpreted cautiously and explicitly as scientifically meaningful in its own right (per Section 31): discrete residue-patch clustering and continuum Poisson-Boltzmann electrostatics are different physical abstractions, and a weak correlation here is not treated as a methodological failure of either.

## X. High-coverage sensitivity

`results/tables/stage8_high_coverage_sensitivity.csv`: all four primary PSH/PPC/PNC group differences preserve direction and approximate magnitude in the ≥0.90-coverage subset (1652/1962) versus the full cohort — conclusions do not materially change.

## Y. Noncanonical audit

`results/audits/stage8_noncanonical_residues.csv`: 115 noncanonical residue occurrences (PCA, UNK, MSE, MLY, TYS, SLL, CSX, TPO) across the cohort; **all** correctly carry `kyte_doolittle_norm = NaN` (never invented), `charge_class = 0.0` (neutral default per the source table's own implicit omission).

## Z. Manual validation result

`results/audits/stage8_manual_validation.csv`: 10 VHH + 10 VH (deterministic, sorted by INSTANCE) + high/low PSH outliers + high PPC/PNC outliers + strong positive/negative APBS cases, including the one known real APBS failure (`pdb_00001zmy-A`, correctly shown as `NaN`, not zero). All patch scores and APBS values fall within plausible cohort ranges; no anomalies found.

## AA. Full repository pytest result

**239/239 PASSED** (up from the 205-test pre-Stage-8 baseline; 34 new Stage-8 tests, zero regressions in Stages 2-7 tests).

## AB. Generated figures

`results/figures/stage8/`: fig1 (hydrophobic patch burden, MAIN_CANDIDATE), fig2 (charge patches, MAIN_OR_SUPPLEMENTARY), fig3 (patch topology, SUPPLEMENTARY, explicitly labeled project-specific), fig4 (APBS potential, MAIN_OR_SUPPLEMENTARY_ELECTROSTATICS). Manifest and README included.

## AC. 3D viewer paths

`results/structures/interactive/stage8_vhh_surface_pdb_00008dqu_pdb_00008dqu-A.html` (VHH representative, blue) and `stage8_vh_surface_pdb_00004ij3_pdb_00004ij3-C-B.html` (VH representative, purple) — deterministic selection near joint cohort medians of `former_interface_core_PSH` and `interface_potential_absolute_mean` (`stage8_representative_selection.csv`). Three separate, non-merged coloring modes (hydrophobic / discrete charge / APBS-sampled-at-atoms) per structure; APBS mode explicitly labeled "sampled at exposed interface atoms," never "surface potential."

## AD. Confirmation: Stage 2-7 scientific outputs unchanged

Verified directly: Stage-4 VHH representative count = 533 (unchanged); Stage-6 combined table N = 533, `cdr3_theta`/`cdr3_phi` still 100% null (unchanged, unused). No Stage 2-7 file was modified this session.

## AE. Confirmation: Gordon theta/phi remain unresolved and unused

Confirmed by direct test (`test_gordon_theta_phi_not_invoked`) across all Stage-8 output tables — no `cdr3_theta`/`cdr3_phi` column present anywhere in Stage-8 outputs; Stage-6's own frozen values remain 100% null, untouched.

## AF. Limitations (Section 58, restated)

1. Patch and electrostatic metrics depend on the chosen (crystallographic) structural conformation.
2. Bound-state crystal structures may bias surface organization and exposure.
3. Crystal packing may affect exposure and electrostatics.
4. TAP/TNP patch descriptors are empirical spatial descriptors, not free energies.
5. Adapted former-interface metrics are project-specific applications of the validated patch framework — never re-derived formulas.
6. PPC/PNC are discrete charge-patch descriptors, distinct from continuum electrostatics.
7. APBS is a continuum electrostatic (linearized Poisson-Boltzmann) calculation, not an exact solution.
8. PDB2PQR protonation assumptions (PARSE force field, propka at pH 7.4) affect all APBS results.
9. `VH_WITHOUT_VL_GEOMETRIC` is a same-coordinate geometric counterfactual, never an experimentally observed isolated VH.
10. Stage 8 does not establish experimental solubility, stability, or aggregation resistance.
11. APBS potential is sampled at exposed atom coordinates, not a true molecular-surface mesh integral — labeled honestly throughout as "potential sampled at exposed interface atoms."
12. 51/1962 structures failed APBS due to a genuine PQR fixed-width-format parsing limitation in this AMBER-era APBS build, unrelated to molecular chemistry; tracked explicitly, never zeroed.
13. Standard TAP/TNP author-code output-level cross-validation was blocked by a real third-party code defect (documented), so standard-descriptor validation rests on source-formula fidelity plus the real `psa` binary, not an author-code output match.

---

## L. Component-aware Stage-8 statuses (reassessed, this closure audit) and M. rationale

Applying the **predeclared** primary/secondary hierarchy literally (primary mechanistic components: adapted former-interface patch analysis, project-specific hydrophobic topology, APBS electrostatics; secondary/contextual: standard whole-antibody TAP/TNP reproduction, SFvCSP) — not defaulting to the prior report's PARTIAL, and not forcing PASS either:

**PATCH STANDARD-METRIC STATUS: PARTIAL** *(secondary component)* — formula source-verified (RESOLVED, taxonomy above) and computed; author-code output-level cross-validation blocked by a genuine, documented third-party `AUTHOR SOFTWARE EXECUTION BUG` (not a method-definition gap); SFvCSP unresolved and correctly not computed.

**ADAPTED FORMER-INTERFACE PATCH STATUS: PASS** *(primary component)* — source-verified formula, correctly retargeted to the frozen empirical interface, validated on 20 real structures, full 1962-structure cohort computed with zero failures.

**PROJECT-SPECIFIC TOPOLOGY STATUS: PASS** *(primary component)* — explicitly labeled non-TAP/TNP, primary threshold reused directly from source without contradiction, 3-threshold sensitivity performed, core finding (largest-patch size/dominance) stable across all thresholds; the conceptual distinction from Stage 5 (burden vs spatial organization) is explicit.

**APBS ELECTROSTATICS STATUS: PASS** *(primary component)* — environment verified, smoke test + validation gate + runtime benchmark all passed on real structures including multi-char-chain/insertion-coded cases; full 1962-structure cohort run with exact, non-ambiguous accounting (1907/1907 successes fully sampled, 0 grid-invalid, 0 sampling-failures among successes, 55 explicit non-success failures never zeroed); interpretation corrected to report positive and negative potential components separately (Section U); honest atom-sampling labeling maintained and re-confirmed (Section D); PPC/PNC and APBS never merged (Section W, weak correlation reported as scientifically meaningful, not as a failure).

**Rationale for the reassessment:** all three PRIMARY components are valid, complete, and tested (verified fresh against frozen outputs in this closure audit, not merely re-asserted). The two remaining limitations — standard TAP/TNP author-code reproduction and SFvCSP — are both explicitly SECONDARY per the original Stage-8 specification's own framing, and both are genuine, honestly-documented limitations (a third-party execution bug and a genuinely-absent external method definition) rather than gaps in this project's own primary analysis. Per the literal predeclared rule ("if all primary components are valid, complete, and tested, and only secondary SFvCSP / standard author-code reproduction limitations remain, those limitations should not automatically force the overall stage to PARTIAL"), retaining PARTIAL at the overall level would overstate the impact of two secondary, non-blocking limitations on a fully-completed primary mechanistic analysis.

## OVERALL STAGE 8 STATUS: PASS

All three primary Stage-8 components (adapted former-interface PSH/PPC/PNC, project-specific patch topology, APBS electrostatics) are complete, source-faithful where applicable, validated on real structures, run across the full 533+1429 cohort with exact and honest failure accounting, and covered by 34 passing tests. The two secondary limitations (standard TAP/TNP author-code reproduction blocked by a documented third-party bug; SFvCSP genuinely unavailable and correctly left NA/unresolved rather than fabricated) are reported in full as limitations (Section AF) but do not, per the predeclared hierarchy, downgrade the overall stage outcome.

## VISUALIZATION DECISION: REQUIRED (fulfilled — 4 figures + 2 interactive 3D viewers with 3 separate coloring modes each)

STOP before Stage 9. No Stage-9+ metrics (antigen engagement, new paratope analyses, aggregation propensity, AggreScan3D, developability classifiers, formal mixed-effects inference, MD) were computed in this stage.
