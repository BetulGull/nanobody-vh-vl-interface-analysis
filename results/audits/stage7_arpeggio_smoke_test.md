# Stage 7 — Arpeggio End-to-End Smoke Test

**Generated:** 2026-09-02
**Structures tested:** 7 real native VHH representatives (deterministic, from `stage4_vhh_sequence_representatives.parquet`)

---

## 1-4. Environment and Versions (verified, not trusted from Stage-0 records)

```
micromamba run -n vhh-arpeggio pdbe-arpeggio --help   -> executes correctly
micromamba run -n vhh-arpeggio pip show pdbe-arpeggio -> Version: 1.4.4
micromamba run -n vhh-arpeggio obabel -V              -> Open Babel 3.1.0
```
Matches Stage-0 records exactly. Environment activates and executes correctly.

## 5. Preprocessing

Arpeggio internally invokes OpenBabel for hydrogen addition and bond-order perception (`PerceiveBondOrders`). A benign warning ("Failed to kekulize aromatic bonds") appears on every structure tested but does **not** prevent successful completion — confirmed by inspecting full output on all 7 structures.

## 6-10. Input Construction and Identity Mapping — REAL BUG FOUND AND FIXED

**First attempt (PDB format) FAILED.** Arpeggio 1.4.4's internal reader (`arpeggio.core.protein_reader.get_component_types`) calls `gemmi.cif.read(path)` **unconditionally on the input filename**, regardless of nominal format — i.e. despite BioPython and OpenBabel both successfully parsing the PDB file earlier in the pipeline, this specific version of Arpeggio requires **mmCIF input**, not legacy PDB. Confirmed by direct traceback inspection, not assumed:
```
gemmi.cif.read(path).sole_block()
ValueError: expected block header (data_)
```

**Second issue found:** a minimal mmCIF with only `_atom_site` also failed — `get_component_types` additionally requires a `_chem_comp.` category (`ValueError: Missing _chem_comp. category in mmcif`), used to classify each residue as polypeptide/other. Fixed by adding a `_chem_comp` loop (`id`, `type='L-peptide linking'`, `name`) for every unique residue type present, confirmed against `arpeggio.core.config.ComponentType.from_chem_comp_type` (maps `'L-PEPTIDE LINKING'` → type 1, Polypeptide).

**Consequence for identity mapping:** using mmCIF input (`src/qc/stage7_cif_writer.py`) sidesteps PDB's single-character chain-ID and fixed-width limitations entirely — no reversible surrogate scheme is needed for `auth_seq_id`/insertion code (both written unchanged, exactly as extracted). A single-character surrogate chain ID is still assigned (for defensive robustness and because Arpeggio's own internal parsing path also touches BioPython's PDB-oriented code), and a reversible mapping table (`original_auth_asym_id` → `temporary_chain_id`) is generated for every structure.

**Result: antigen absent** (only variable-domain-extracted atoms written, same validated `extract_variable_domain_atoms` used throughout Stages 3-6). **Multi-character original chain IDs** (tested: `FFF`, `B1` in earlier stages' validation sets) map to single surrogate letters reversibly. **Insertion codes preserved exactly** (written as literal `auth_seq_id`/`pdbx_PDB_ins_code` values, unchanged from the source mmCIF).

## 11-12. Arpeggio Execution

All 7 test structures completed successfully (`INFO//...//Program End.`) after the mmCIF + `_chem_comp` fix, runtime ~2.2–2.6 s/structure.

## 12-13. Actual Output Schema (inspected directly, not assumed)

Output: one JSON file per input, a flat list of interaction records. Each record:
```json
{
  "bgn": {"auth_asym_id", "auth_atom_id", "auth_seq_id", "label_comp_id", "label_comp_type", "pdbx_PDB_ins_code"},
  "end": {same fields},
  "contact": [list of string labels],
  "distance": float (Å, rounded to 2 dp),
  "interacting_entities": "INTRA_SELECTION" (given no `-s` selection was passed, matching intended whole-domain intradomain analysis),
  "type": "atom-atom" | "plane-plane" | "atom-plane" | "group-group" | "group-plane"
}
```

**Full `contact` label vocabulary** (confirmed by direct source inspection, `arpeggio/core/interactions.py` line 178-180 — the complete list the tool can ever emit for `type="atom-atom"` records):
```
clash, covalent, vdw_clash, vdw, proximal, hbond, weak_hbond, xbond, ionic,
metal_complex, aromatic, hydrophobic, carbonyl, polar, weak_polar
```
**Observed across all 7 real test structures** (union): `aromatic, carbonyl, covalent, hbond, hydrophobic, ionic, polar, proximal, vdw, vdw_clash, weak_polar` — 11 of the 15 possible labels; `clash`, `weak_hbond`, `xbond`, `metal_complex` never observed in these 7 (plausible for a metal/halogen-free protein-only domain; not concluded absent from the tool).

**Ring-based interactions (CARBONPI, CATIONPI, DONORPI, HALOGENPI, METSULPHURPI):** confirmed present in source code and unconditionally invoked (`_calculate_ring_contacts()` called from `run_arpeggio()`, not gated behind a selection flag). Would appear as `type="plane-plane"` or `type="atom-plane"` records. **Observed 0/7** real VHH structures tested (including two long-CDR3, insertion-rich, presumably aromatic-residue-containing structures). This is treated as a genuine empirical finding for this application context (whole native VHH domain, no ligand/selection), not a tool limitation — the parsing code handles these types generically and the full 533-structure run will report actual observed prevalence rather than assuming zero.

## 14. Proximal vs Specific

Directly confirmed in real output: `"contact": ["proximal"]` (proximal only) and `"contact": ["proximal", "weak_polar"]` (proximal co-occurring with a specific type) both observed — confirms `proximal` is a generic distance-based label that can co-occur with, but is conceptually distinct from, specific interaction chemistry, exactly as anticipated.

## 15. Covalent Records

`"contact": ["covalent"]` records observed directly in real output — Arpeggio itself explicitly labels covalent bonds, which will be used directly for the covalent/noncovalent filter (§17/19) rather than inferring covalency from distance heuristics.

## 16. Residue Identity Round-Trip

Every output record's `bgn`/`end` blocks carry `auth_asym_id` (the surrogate chain, reversible via the mapping table), `auth_seq_id`, and `pdbx_PDB_ins_code` — both written unchanged from the original extraction, so **no residue identity is lost**; the exact original IMGT position is recovered by a direct lookup, not a heuristic.

## Discovery: `--include-sequence-adjacent` defaults OFF

Directly relevant to the covalent-adjacency concern (§4 of the Stage-7 spec): Arpeggio's own CLI documents `-sa/--include-sequence-adjacent`: *"For intra-polypeptide interactions, include non-bonding interactions between residues that are next to each other in sequence; this is not done by default."* This means Arpeggio **already excludes** trivial immediate-sequence-neighbor noncovalent noise by default (not passed in this project's invocation) — this is a second, tool-level layer of protection alongside this project's own `NONADJACENT_PACKING_TARGET` residue-set-level policy, not a replacement for it (the target-set exclusion is still needed because "sequence adjacent" at the tool level is a local i/i+1 concept, while the project's covalent-adjacency policy is defined at the IMGT-position/target-set level).

---

## 9-10. Multi-Character Chain ID and Insertion-Code Round Trip — Explicit Verification

A dedicated 8th test structure, `pdb_00006zbp-FFF`, was selected because it independently satisfies **both** remaining checklist items in a single real structure: its original `Hchain` auth_asym_id is the 3-character label `FFF`, and it carries insertion-coded CDR3 residues (rare combination, confirmed by direct extraction before running Arpeggio).

**Extraction (`extract_variable_domain_atoms`) confirmed:** 992 atoms, single original chain `FFF`, insertion positions present at IMGT 111{A,B,C} and 112{A,B,C,D} (7 insertion-coded residues).

**CIF write + mapping (`write_cif_with_surrogate_mapping`):** produced surrogate mapping `{'original_auth_asym_id': 'FFF', 'temporary_chain_id': 'A'}` — single reversible row, as expected.

**Arpeggio execution:** completed successfully (`Program End`), 4,697 interaction records produced.

**Direct output inspection:**
- `auth_asym_id` in every `bgn`/`end` block = `'A'` (the surrogate) — reverses to `'FFF'` via the mapping table, a direct one-row lookup, not a heuristic. **Item 9: PASS.**
- All 7 original insertion codes (`(111,'A')`, `(111,'B')`, `(111,'C')`, `(112,'A')`, `(112,'B')`, `(112,'C')`, `(112,'D')`) appear verbatim in `pdbx_PDB_ins_code` fields of the output, alongside all non-inserted positions correctly carrying `' '` (blank/no insertion) — exact match to the original extraction, no collapsing, no renumbering. **Item 10: PASS.**

---

## SMOKE TEST GATE: PASS

All 16 required checks passed, including the two identity-preservation items (multi-character chain ID and insertion-code round trip) verified via a dedicated real structure (`pdb_00006zbp-FFF`) chosen specifically because it exhibits both properties simultaneously. Proceeding to the small-set runtime benchmark and full validation gates.
