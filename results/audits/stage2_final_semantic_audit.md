# Stage 2 — Final Semantic Audit

**Generated:** 2026-09-01
**Scope:** Verification-only pass over the closure-audit-fixed Stage-2 pipeline. No implementation change was ultimately required; no full-cohort rerun performed.

---

## 1. Expected-Residue → IMGT Mapping Source

`_pdbx_poly_seq_scheme` confirmed absent (pre-checked, not rediscovered). No independent scheme table exists in these files for the IMGT/auth number of a completely unresolved residue.

**What was validated** (`results/audits/stage2_expected_residue_imgt_mapping_validation.csv`, 5 deterministic cases: legitimate IMGT gap, fully complete CDR, insertion-coded CDR3, true missing-coordinate residue, multiple missing residues in one CDR):
- the true polymer (`_entity_poly_seq`) contains additional real residues beyond the coordinate block, exactly when that's actually true
- those residues lie between correctly-identified observed anchors, with count and order matching the polymer's `label_seq_id` sequence
- legitimate IMGT numbering gaps are never counted as missing (0 false positives across all 5 cases)
- insertion-coded residues (e.g. 112, 112A, 112B, 112C, 112D in the 9-missing-residue example) remain distinct and correctly ordered

**Result: 5/5 cases pass every provable check.**

**Explicit limitation (not overclaimed):** the exact IMGT/auth number a completely unresolved residue "would have had" cannot be independently proven without a scheme table this dataset doesn't have. The method never fabricates one — `get_missing_cdr_positions()` reports `auth_seq_id=None` for such residues, keyed only by `label_seq_id` and monomer identity.

**Granular missingness counts** (`results/audits/stage2_cdr_missingness_breakdown.csv` + `_summary.csv`, replacing the prior vague "163 flags"):

| Metric | Count |
|---|---|
| Total missing real CDR residues | **678** |
| Affected CDR loops (INSTANCE, CDR pairs) | **189** |
| Affected antibody chains (INSTANCE, chain pairs) | **179** |
| Affected INSTANCE rows | **173** |

Per-CDR breakdown (residues missing / loops affected):

| CDR | Residues missing | Loops affected |
|---|---|---|
| H1 | 217 | 72 |
| H2 | 40 | 16 |
| H3 | 192 | 50 |
| L1 | 121 | 30 |
| L2 | 82 | 11 |
| L3 | 26 | 10 |

---

## 2. Antibody-Chain Antigen Masking Semantics

Every chain suppressed by the multi-instance mask (576 suppressed-chain instances across the 1244 multi-INSTANCE PDBs / 3235 INSTANCEs) was classified using: this INSTANCE's own SAbDab2 `antigen_chain`/`antigen_type`/`antigen_name` fields, the owning INSTANCE's `SABDAB_ID` relationship, and the owning INSTANCE's own `antigen_chain` field (mutual cross-reference) — full detail in `results/audits/stage2_antibody_antigen_mask_semantics.csv`.

**First pass finding (important):** a naive cross-check ("does this INSTANCE's own `antigen_chain` field name the suppressed chain?") returned **576/576 "legitimate"** — which was a red flag, not a result. Manual inspection of the underlying examples showed SAbDab2's own `antigen_chain` field is itself unreliable for multi-INSTANCE PDBs, in exactly the way the closure-audit fix anticipated:
- **Apo antibody-dimer self-reference:** e.g. `pdb_00003qxu` chains B and D — antigen_name "Anti-Methotrexate CDR1-3 Graft" describes the antibody *itself*, mutually mislabeling its own crystallographic dimer partner as "antigen."
- **Multi-antibody complexes against a shared real target:** e.g. `pdb_00006fe4`, antigen_name "Shiga-like toxin 2 subunit B|Nb113" — SAbDab2 lumps a co-bound sibling nanobody (Nb113) into the antigen_name/antigen_chain fields alongside the genuine target.
- **Named historical antibody as ternary-complex bystander:** `pdb_00009cds`, antigen_name "Outer surface protein A|Hybridoma Antibody LA2" — a third, specifically-named mAb (LA2) bound to the same real antigen (OspA) at a different epitope, not a deliberate antibody-as-target relationship for the INSTANCE being scored.

**Refined criteria** (documented in `stage2_antibody_antigen_masking_audit.py`) require a strong anti-idiotype/antibody-binding-protein signal (`ANTI-IDIOTYPE`, `IDIOTOPE`, `PROTEIN A/G/L`, `FC RECEPTOR`, `RHEUMATOID FACTOR`) **with no sibling-antibody contamination** (`NANOBODY`, `FAB`, `SCFV`, `VHH`, `HYBRIDOMA`, `MONOCLONAL ANTIBODY`, named light/heavy chain fragments, etc.) in the antigen_name text before accepting a chain as legitimately antigen-like, applied symmetrically to both this INSTANCE's own field and the owning INSTANCE's mutual cross-reference.

**Final classification breakdown:**

| Classification | Count |
|---|---|
| SAME_ANTIBODY_COMPONENT (same `SABDAB_ID` — literally the same antibody sequence/entity) | 330 |
| SIBLING_CRYSTAL_COPY_OR_PARALLEL_ANTIBODY_INSTANCE (ordinary co-deposited antibody, own separate `antigen_chain`) | 246 |
| POTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN | **0** |
| UNRESOLVED | 0 |

**Result: zero legitimate antibody-derived-antigen cases found among all 576 suppressed-chain instances.** Every suppressed chain is either (a) a duplicate/component of the very same antibody sequence (`SABDAB_ID` match) or (b) an ordinary separate antibody targeting its own real antigen in the same crystal. Spot-checks of one example from each category confirm the classification logic is sound (`pdb_00003ezj-B`↔F: same `SABDAB_ID`, correctly `SAME_ANTIBODY_COMPONENT`; `pdb_00005ocl-A`↔B: different `SABDAB_ID`, distinct own `antigen_chain`, correctly `SIBLING_CRYSTAL_COPY`).

**Conclusion: the current global multi-instance antibody-chain mask is empirically justified for this dataset. No refinement to the masking logic is required.**

---

## 3. Current Final Counts — Confirmed Unchanged

No implementation change resulted from either check, so no rerun was performed (per instruction). Direct verification against the saved tier files:

| Cohort | Count |
|---|---|
| Strict VHH structural cohort | **1217** |
| Paired primary structural cohort | **2879** |
| VHH secondary/sensitivity structural cohort | **126** |

Full pytest suite: 64/64 passed (unchanged, no code modified in the pipeline itself — only new audit scripts added).

---

## STAGE 2 FINAL SEMANTIC AUDIT: PASS
