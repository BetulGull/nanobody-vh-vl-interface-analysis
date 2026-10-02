# Stage 9 — Antigen Engagement, Framework Contribution, and Functional Accessibility of the Former VH-VL Surface

## A/B. Starting cohorts

**Starting VHH N = 533**, **starting conventional N = 1429** — the frozen Stage-3/4 primary nonredundant representative cohorts, unchanged.

## C/D/E. Protein-antigen-engaged cohort derivation and exclusion funnel

Per the amendment, the Stage-2 antigen/Gordon-contact gate was traced back to its generating code (`stage9_antigen_engagement_method_definition.md`) and verified to faithfully implement every required Stage-9 sub-component (legitimate protein-antigen identification via mmCIF entity type, sibling-antibody and multi-instance-PDB masking, >=1 CDR-residue <=4.5Å heavy-atom contact, multi-antigen-chain handling). **Reused rather than recomputed.**

**Finding, verified directly against the frozen files (not assumed):** both frozen primary representative cohorts are already 100% antigen-engaged:

```
stage4_vhh_sequence_representatives.parquet (N=533): qc_gordon_cdr_antigen_contact == True: 533/533
stage3_paired_representatives.parquet (N=1429, primary): qc_gordon_cdr_antigen_contact == True: 1429/1429
```

- **Protein-antigen-engaged VHH N = 533** (100% of starting cohort)
- **Protein-antigen-engaged conventional N = 1429** (100% of starting cohort)

**Exclusion funnel** (full Stage-2 candidate pool, for context — not the representative-level funnel, since no further reduction occurs at the representative level): VHH 1380 core-structure-pass → 1342 Gordon-complex-pass (38 excluded upstream); paired 2997 core-structure-pass → 2879 Gordon-complex-pass (118 excluded upstream). The downstream Stage-3/4 redundancy-reduction/representative-selection steps drew exclusively from the Gordon-passed pool, so the final 533/1429 representatives require no additional Stage-9 filtering.

## F. Antigen-identification validation

`AntigenQuality.identify_antigen_chains` traced and verified: mmCIF-entity-type-based (`polypeptide(L)`/`polypeptide(D)` only), explicit exclusion of water/ion/non-polymer/nucleic-acid/glycan, sibling-antibody and multi-instance-PDB chain masking via `build_pdb_antibody_chain_map` (all SAbDab2 instances per PDB, not just the current row).

## G. Contact-definition validation

20-structure (10 VHH + 10 VH) cross-validation against **independently-computed Stage-2 ground truth** (`contacting_cdrs`, from a completely separate code path, `contact_search.py`): **20/20 exact match** on which CDRs (H1/H2/H3/L1/L2/L3) register antigen contact (`results/audits/stage9_antigen_contact_validation.csv`). This is a stronger validation than a generic spot-check — it confirms the new Stage-9 residue-pair-level contact machinery reproduces the already-frozen, independently-derived Stage-2 result exactly.

**A real data-cleanliness bug was found and fixed during this validation work:** `AntigenQuality.get_antigen_atoms` returns all atoms sharing an antigen chain's author chain ID, including crystallographic waters (`HOH`, `group_PDB=='HETATM'`) that happen to share that chain ID — a genuine, common mmCIF pattern. This was initially silently included in the antigen atom set (confirmed: 34,095/152,952 = 22.3% of initial contact-pair rows were water, not real protein-antigen contacts), discovered via a real Arpeggio smoke test surfacing an unexpected `HOH` residue in the antigen selection. **Fixed** by filtering to `group_PDB=='ATOM'` at the atom-record level (never touching the frozen Stage-2 entity-level chain filter), and the **full 1962-structure primary cohort was rerun** with the fix. The corrected long table (118,614 rows) contains zero water contacts. The Stage-2 cross-validation was re-run post-fix and still shows 20/20 exact match — water inclusion never changed the boolean CDR-contact classification in this sample, confirming the primary conclusions are robust to this correction, but the fix materially improves data cleanliness for all downstream residue/atom-level counts.

## H. Antigen composition audit

`results/audits/stage9_antigen_composition_audit.csv`: no single antigen dominates either cohort (VHH top antigen "Spike protein S1" = 6.9% of cohort; VH top antigen "Spike protein S1" = 7.3%) — conclusions are not driven by one repeated antigen family.

## I/J. Total paratope size

`results/tables/stage9_vhh_vs_fv_total_paratope.csv`: VHH total unique antibody contact residues, unique antigen contact residues, and unique residue pairs, compared against the FULL conventional Fv (VH+VL) paratope — the correct total-paratope comparison per Section 8's three-tier framework.

## K/L. Framework contribution — VHH vs conventional VH (VH-only, not Fv)

`results/tables/stage9_vhh_vs_vh_framework_engagement.csv`:

- **VHH framework_contact_fraction: mean 0.209** vs **conventional VH framework_contact_fraction: mean 0.163** — VHH shows greater framework participation in the observed antigen interface than the conventional VH domain alone.
- Determined empirically which framework regions participate (`stage9_framework_antigen_contact_frequency.csv`, 77 framework positions show any contact): the largest VHH-vs-VH differences concentrate in FR2 (positions 52, 42, 50, 49 — the hallmark tetrad region) and FR4 (position 118), but FR3 and FR1 also contribute meaningfully (33 and 22 of the 77 contacting positions respectively) — **FR2 does not exclusively dominate**, confirmed empirically rather than assumed.

## M/N. CDR3 contribution

`results/tables/stage9_cdr3_antigen_contribution.csv`: **VHH cdr3_contact_fraction: mean 0.470** vs **conventional VH: mean 0.403** — VHH antigen engagement is somewhat more CDR3-dominated than conventional VH, consistent with (but not proof of) the well-known extended-CDR3 VHH paradigm. No causal claim made.

## O/P. Former-interface antigen engagement (VHH) — full 28-position and 15-position non-CDR3

`results/tables/stage9_former_interface_antigen_engagement.csv`:

- **Full 28-position empirical interface: VHH mean engagement fraction 0.357** — over a third of the former VH-VL interface positions contact antigen on average.
- **15-position non-CDR3 shieldable interface: VHH mean engagement fraction 0.186** — even excluding all CDR3 positions, nearly a fifth of the purely-framework former-interface positions contact antigen.
- Both kept strictly separate throughout (never conflating CDR3 antigen contacts with framework former-interface engagement, per Section 16).

## Q. Conventional VH former-interface antigen engagement

Same table: **VH full 28-position: mean 0.224**; **VH 15-position non-CDR3: mean 0.089**. VHH engagement is higher at the former interface for both target definitions than the corresponding conventional VH positions (which, in conventional antibodies, are normally buried against VL).

## R. Tetrad antigen engagement

`results/tables/stage9_tetrad_antigen_engagement.csv` (IMGT 42, 49, 50, 52 — annotation-only, never redefining the former interface): dramatic VHH-vs-VH differences at every tetrad position (e.g., position 52: VHH 34.0% vs VH 9.2%; position 42: VHH 22.9% vs VH 0.3%) — these hallmark positions are substantially functionally engaged with antigen in VHH, far more than their structural counterparts in conventional VH, where they are normally packed against VL rather than solvent/antigen-facing.

## S. Shielding vs antigen-engagement relationship (VHH only) — central Stage-9 finding

`results/tables/stage9_shielding_vs_antigen_engagement.csv` (N=533, all Pearson/Spearman, descriptive/exploratory only):

| relationship | Pearson r |
|---|---|
| shielding_fraction vs full-interface engagement | -0.29 |
| shielding_fraction vs 15-position non-CDR3 engagement | -0.44 |
| shielding_fraction vs framework_contact_fraction | -0.40 |
| shielding_fraction vs cdr3_contact_fraction | +0.42 |

**A consistent, moderate COMPETITION-like pattern**: more CDR3 shielding is associated with less former-interface/framework antigen engagement and more CDR3-dominated engagement. This is reported descriptively — as a measured relationship, not an assumed "cost of shielding" — and the correlations are moderate (not overwhelming), leaving room for both competition and coexistence within the cohort (manual validation, Section Z, confirms individual cases of both high-shielding-high-engagement and high-shielding-low-engagement genuinely exist).

## T. Packing vs antigen-engagement relationship

`results/tables/stage9_packing_patch_vs_antigen_engagement.csv`: Stage-7 `nonadj_BSA_pair_total` (CDR3-target packing) shows the same competition-like pattern (r=-0.33 full interface, r=-0.48 shieldable), extending the shielding finding to native CDR3-target packing.

## U. Geometry vs antigen-engagement relationships

`results/tables/stage9_geometry_vs_antigen_engagement.csv`: CDR3 sequence length (r=-0.24), Gordon compactness (r=-0.25), and Zhang FR2_CDR3_Dis (r=+0.26) all show modest relationships with full-interface engagement, consistent in direction with the shielding finding (more compact/closer CDR3 geometry associated with less former-interface engagement). Gordon ρ shows no meaningful relationship (r=0.06, p=0.19) — not every CDR3 geometry descriptor relates to engagement. Gordon θ/φ not used (remain unresolved, per Stage 6).

## V. Stage-8 patch/electrostatic descriptors vs antigen engagement

`results/tables/stage9_packing_patch_vs_antigen_engagement.csv`: Stage-8 `former_interface_core_PSH` shows the same modest inverse relationship (r=-0.28) with shieldable engagement; PPC/PNC show negligible relationships (|r|<0.09) — the competition-like pattern is clearest for hydrophobic packing/patch descriptors, not for discrete charge patches.

## W. Threshold sensitivity

`results/tables/stage9_contact_threshold_sensitivity.csv`: monotonic increase in contact-pair and unique-residue counts from 4.0→4.5→5.0 Å for both groups, as expected; the 4.5 Å primary threshold was not chosen after inspecting group differences (predeclared, source-verified from the Gordon criterion already used in frozen Stage-2 code).

## X. High-coverage sensitivity

`results/tables/stage9_high_coverage_sensitivity.csv` (394/533 VHH structures with >=0.90 primary-interface position coverage): all four primary engagement fractions preserve direction versus the full cohort — conclusions do not materially change.

## Y. Arpeggio antigen-interface status — SECONDARY LAYER, COMPLETED

Per the amendment's ordering, this layer was attempted only after all primary geometric analyses were complete and validated. Treated as a genuinely new inter-partner use case (Stage 7's intradomain validation was not assumed to transfer):

1. **CLI/API inspection**: `-s SELECTION` accepts `/<chain_id>/<res_num>[<ins_code>]/<atom_name>` syntax; `/A//` selects a whole chain as the "ligand," producing an `interacting_entities` field with three real values — `INTRA_SELECTION` (within the selected chain, i.e. within antibody), `INTRA_NON_SELECTION` (within the non-selected chains, i.e. within antigen), and **`INTER`** (the antibody-antigen cross-interactions of interest) — discovered by direct inspection, not assumed from Stage 7 (which never used `-s`, since Stage 7 was single-partner).
2. **One real VHH-antigen smoke test**: initially failed twice with real, fixable bugs (see below), succeeded on the third attempt (5,688 records, 524 INTER).
3. **One real Fv-antigen smoke test**: succeeded using the same fixed pipeline.
4. **Inter-partner schema discovered directly**: same `bgn`/`end`/`contact`/`distance`/`type` record shape as Stage 7, with the new `interacting_entities` values as the key addition.
5. **Antibody<->antigen residue identity round-trip validated**: surrogate-chain reverse-mapping applied to every INTER record's antigen-side chain ID; antibody IMGT position, insertion code, and region classification carried through unchanged.
6. **One genuine multi-instance PDB verified**: sibling-instance contamination is structurally impossible by construction — the combined CIF for a given instance contains only that instance's own antibody chain(s) and its own Stage-2-assigned antigen chains, never any other SAbDab2 instance's atoms.
7. **5 VHH + 5 Fv validation gate**: 10/10 succeeded after fixes (`results/audits/stage9_arpeggio_antigen_validation.csv`), 341-1364 INTER records per structure, runtime 2.8-4.2s.
8. **Runtime benchmark**: ~3.5s/structure mean (10 real structures spanning antigen sizes 1148-3909 atoms) — full 1962-structure run projected at ~2 hours.
9. **Production decision**: full cohort judged practical (comparable cost/structure to Stage 7/8's successful full runs).
10. **Full antigen-engaged cohort run**: **1962/1962 `ARPEGGIO_SUCCESS`**, zero failures.

**Two real, fixable bugs found and fixed during this work (documented, not silently patched around):**
- **Bug 1**: a duplicate-atom crash (`Atom CB defined twice`) traced to unfiltered altloc conformers in antigen atom extraction — fixed by applying the same highest-occupancy/alphabetical-tiebreak altloc-selection policy already used for antibody atoms throughout this project (`select_altloc` pattern). This fix also corrected the **primary** geometric contact pipeline (see Section G) since both use the same antigen-atom extraction function — though the primary contact/distance results were unaffected in practice (minimum-distance search is robust to near-identical duplicate atom copies).
- **Bug 2**: an output-path duplication bug in this project's own subprocess invocation (passing an absolute `-o` path while also setting `cwd` to the same directory caused Arpeggio to write output to a spuriously nested directory tree) — fixed by passing `-o .` when `cwd` is already set to the target directory.

**Results** (`results/tables/stage9_arpeggio_specific_interaction_summary.csv`, 1,476,750 total antibody-antigen INTER atom-pair records): 83.7% proximal-only, 16.3% specific interactions (hbond/ionic/aromatic/hydrophobic/polar/vdw etc.), negligible covalent (2.7e-6, as expected for a non-covalent antibody-antigen interface — confirms proximal is correctly distinguished from specific throughout, per the standing Stage-7 convention reused here). Per-structure specific-interaction-atom counts: VHH mean 64.9 (range 2-136), conventional Fv mean 144.8 (range 64-248) — the larger Fv paratope naturally registers more specific-interaction atoms in aggregate, consistent with Section I/J's total-paratope-size finding. Arpeggio counts are explicitly never treated as binding energy (verified by test), consistent with the Stage-7 convention.

**Status: Standard/secondary Arpeggio antigen-interaction layer — COMPLETE (not blocked).**

## Z. Manual validation

`results/audits/stage9_manual_validation.csv` (46 rows: 10 low/10 near-median/10 high VHH framework engagement, 3 high-shield-high-engagement, 3 high-shield-low-engagement VHH cases, 5 low/5 high VH framework engagement): confirms both `high_shield_high_engagement` and `high_shield_low_engagement` VHH structures genuinely exist in the cohort (e.g. `pdb_00001ri8-A`: shielding_fraction=0.31, shieldable_engagement=0.33 vs `pdb_00001jtp-A`: shielding_fraction=0.28, shieldable_engagement=0.00) — confirming Section S's relationship is a real, moderate population-level tendency, not a deterministic rule, and that COEXISTENCE cases are not rare exceptions.

## AA. Full repository pytest result

**270/270 PASSED** (up from the 239-test pre-Stage-9 baseline; 31 new Stage-9 tests, zero regressions in Stages 2-8 tests, re-confirmed after the water-contamination fix and full-cohort rerun).

## AB. Generated figures

`results/figures/stage9/`: fig1 (framework contribution, MAIN_CANDIDATE_FUNCTIONAL_MECHANISM), fig2 (region composition, MAIN_OR_SUPPLEMENTARY_PARATOPE_ORGANIZATION), fig3 (former-interface position map integrating Stage-3 VL burial + VHH/VH antigen-contact frequency, MAIN_CANDIDATE_INTEGRATED_MECHANISM), fig4 (shielding vs engagement, VHH only, title dynamically reflects the observed inverse relationship). Manifest and README included.

## AC. Viewer paths

`results/structures/interactive/stage9_vhh_antigen_pdb_00008pe2_pdb_00008pe2-B.html` and `stage9_fv_antigen_pdb_00002oqj_pdb_00002oqj-B-A.html` — deterministic representatives near joint cohort medians of `framework_contact_fraction` and `shieldable_engagement_fraction`. VHH viewer: blue backbone, cyan CDR3, orange sticks (full former interface), red sticks (framework residues actually contacting antigen), restrained green (native antigen). Fv viewer: purple VH, pink VL, red sticks (VH framework antigen contacts), restrained green antigen.

## AD. Confirmation: Stage 2-8 frozen outputs unchanged

Verified directly: Stage-4 VHH N=533, Stage-3 paired primary N=1429, Stage-8 patch descriptors N=1962 — all unchanged. No Stage 2-8 file was modified this session.

## AE. Confirmation: Gordon theta/phi remain unresolved and unused

Verified by direct test (`test_gordon_theta_phi_not_invoked`) across all Stage-9 output tables — absent everywhere, per Stage 6's frozen PARTIAL status.

## AF. Limitations (Section 56, restated)

1. Stage 9 uses antigen-bound experimental structures — describes observed bound-state paratope geometry, not unbound conformational accessibility.
2. Contact participation is not binding affinity; no KD/kon/koff/ΔG language used anywhere (verified by test).
3. Bound conformations may differ from unbound conformations.
4. Crystal packing may influence CDR3 and interface geometry.
5. Antigen classes are heterogeneous (composition audit performed, no single antigen dominates).
6. Contact cutoffs are geometric conventions (3-threshold sensitivity performed).
7. Arpeggio annotations (where completed) are not interaction free energies.
8. Repeated antigen/epitope states may require formal hierarchical handling in Stage 11 (primary cohort retains frozen nonredundant representatives only).
9. Lack of observed antigen contact does not imply inability to contact another antigen.

---

## STAGE 9 STATUS: PASS

All required conditions verified: protein antigen identification validated by tracing Stage-2 source code and cross-checking 20/20 real structures against independently-computed ground truth; sibling-antibody contamination structurally excluded by construction (Stage 2's multi-instance masking, reused, plus per-instance CIF construction for Arpeggio never including sibling atoms); native antigen-bound states correctly reconstructed (antigen never computationally removed, unlike Stages 5/8); VHH, VH, and VL contributions kept distinct throughout (separate columns, separate test coverage); 4.0/4.5/5.0 Å threshold sensitivity complete and predeclared; full 28-position and 15-position non-CDR3 former-interface definitions kept distinct (never substituting the Stage-7 14-position nonadjacent target); framework/CDR denominators verified correct by test; manual validation passed (46 rows, confirms both competition and coexistence cases genuinely exist); full pytest passes (270/270); Stage 2-8 scientific outputs confirmed unchanged. The secondary Arpeggio antigen-interaction layer, initially expected to be optionally blockable per the status policy, in fact completed successfully on the full 1962-structure cohort after two real, documented, and fixed bugs (one of which also improved the primary geometric layer's data cleanliness, requiring and receiving a full cohort rerun before finalizing).

## VISUALIZATION DECISION: REQUIRED (fulfilled — 4 figures + 2 interactive 3D antigen-interface viewers)

STOP before Stage 10. No Stage-10+ metrics (AggreScan3D, aggregation propensity, developability prediction, formal mixed-effects models, multiple-testing-adjusted final inference, MD) were computed in this stage.
