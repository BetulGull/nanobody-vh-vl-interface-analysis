# Stage 9 — Antigen Engagement Method Definition: Stage-2 Reuse Audit

Per the amendment, this audit traces the frozen Stage-2 antigen/Gordon-contact fields back to their **generating code** (not their names) before any Stage-9 antigen logic is written.

## Traced source: `src/qc/stage2_structure_qc.py` (`process_one_structure`, steps 5-6), calling `src/qc/antigen_qc.py::AntigenQuality` and `src/qc/contact_search.py::ContactSearch`.

### 1. Legitimate protein-antigen identification — VERIFIED EXACT

`AntigenQuality.identify_antigen_chains` (`antigen_qc.py` lines 43-95): classifies every non-antibody chain in the mmCIF `_entity_poly`/`entity_info` mapping. A chain is antigen **only if** its mmCIF entity type is exactly `polypeptide(L)` or `polypeptide(D)` (`PROTEIN_POLYMER_TYPES`). Explicitly excluded (`EXCLUDE_TYPES`): `water`, `ion`, `non-polymer`, `nucleotide`, `nucleic-acid`, `branched` (glycans). This is entity-mapping-based, never inferred from free-text compound/description fields — satisfies Stage-9 Section 3 exactly (protein/amino-acid polymer only; excludes glycans-alone, nucleic-acid-alone, small molecules, water, ions).

### 2. Sibling-antibody / multi-instance PDB masking — VERIFIED EXACT

`stage2_structure_qc.py::build_pdb_antibody_chain_map` (lines 79-100) builds, for every PDB, the union of **every H/L chain used by every SAbDab2 INSTANCE row in both the VHH and paired candidate pools for that PDB** — not just the current row's own H/L pair. `process_one_structure` (line 233) passes `antibody_chains = {hchain, lchain} | pdb_antibody_chains.get(pdb, set())` into `identify_antigen_chains`, which excludes all of them (line 71-72 of `antigen_qc.py`). This means a co-deposited second antibody copy, or any other antibody instance sharing the same PDB file, is **never** misidentified as antigen — directly satisfying Stage-9 Section 3's "exclude sibling antibody chains" / "unrelated antibody instances" and Section 6's multi-chain/multi-instance handling requirement.

### 3. Gordon >=1 CDR-residue, <=4.5 Å heavy-atom contact criterion — VERIFIED EXACT

`ContactSearch` (`contact_search.py`), docstring quotes the exact criterion: *"A CDR is 'in contact with' protein antigen if at least one CDR residue has at least one atom within <=4.5 Å of at least one antigen atom."* Instantiated with `threshold_angstrom=4.5` (`stage2_structure_qc.py` line 266). CDR atoms are extracted per-CDR (H1/H2/H3, +L1/L2/L3 if paired) via `CDRQuality.get_cdr_atom_coords`, restricted to CDR-region residues only (not whole-chain) — so a contact registers only via genuine CDR-region atoms, exactly matching the "at least one CDR residue contacting antigen" requirement. Distance computed via a validated dual implementation (`brute_force_contacts` reference + `optimized_contacts` cKDTree, both documented to return identical results) — heavy-atom (all non-hydrogen structural atoms, since mmCIF `_atom_site` records used throughout this project never include hydrogens at the deposited-structure stage) minimum pairwise distance, not CA-only or centroid-based.

### 4. Multiple antigen-chain handling — VERIFIED

`antigen_chains` (list) and `n_antigen_chains` (count) are retained per-instance; `contacting_antigen_chains` records exactly which of the candidate antigen chains actually registered a contact (line 274) — chains present but not contacting are distinguishable from chains that do contact, satisfying Stage-9 Section 6's "preserve antigen_chain_count and exact identities."

## Verdict: the frozen Stage-2 gate FAITHFULLY implements the Stage-9 Gordon-inspired antigen-engagement gate

All five required sub-components (legitimate protein-antigen ID, >=1 CDR-residue contact, 4.5 Å heavy-atom criterion, multi-antigen-chain handling, sibling-antibody/multi-instance masking) are verified present and correctly implemented in the frozen Stage-2 code and its outputs. **Per the amendment: these frozen results are REUSED for Stage-9 cohort derivation rather than recomputed.**

## Direct verification against the frozen primary representative files (not assumed from column names)

```
stage4_vhh_sequence_representatives.parquet (N=533):
  qc_has_candidate_protein_antigen == True : 533/533 (100%)
  qc_gordon_cdr_antigen_contact    == True : 533/533 (100%)

stage3_paired_representatives.parquet (is_primary_representative==True, N=1429):
  qc_has_candidate_protein_antigen == True : 1429/1429 (100%)
  qc_gordon_cdr_antigen_contact    == True : 1429/1429 (100%)
  qc_stage2_gordon_complex_pass    == True : 1429/1429 (100%)
```

**Finding, reported immediately per instruction (not assumed to be close to 533/1429 in advance, but verified to be EXACTLY 533/1429):** the full frozen Stage-3/4 primary nonredundant representative cohorts are **already 100% antigen-engaged** under the Gordon criterion. This is explained by the Stage-2 funnel (`stage2_cohort_funnel.csv`): Gordon-complex-pass was a real, non-trivial filter at the full candidate-pool level (VHH: 1380 core-pass -> 1342 Gordon-pass, i.e. 38 excluded; paired: 2997 core-pass -> 2879 Gordon-pass, i.e. 118 excluded), and the downstream redundancy-reduction/representative-selection steps (Stages 3-4) that produced the final 533/1429 representatives drew exclusively from structures that had already passed this gate — so no further reduction occurs at the representative level.

**Consequence for Stage 9:** `VHH_ANTIGEN_ENGAGED_PRIMARY = the full frozen 533-VHH cohort` and `CONVENTIONAL_ANTIBODY_ANTIGEN_ENGAGED_PRIMARY = the full frozen 1429-VH cohort`. No new antigen-engagement filtering step reduces these cohorts further — Stage 9's "analysis subset" is, empirically, identical to the Stage-3/4 frozen cohorts. This is reported as a genuine empirical finding of this audit, not assumed in advance.
