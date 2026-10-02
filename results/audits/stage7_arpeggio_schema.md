# Stage 7 — Arpeggio Output Schema (real, source-confirmed)

Derived from direct inspection of real Arpeggio 1.4.4 JSON output across the full 533-structure cohort (`data/interim/stage7/cohort_json/`), cross-checked against `arpeggio/core/interactions.py` and `arpeggio/core/config.py` source. See `results/audits/stage7_arpeggio_smoke_test.md` for the initial 7/8-structure derivation; this document reflects the full-cohort confirmation.

## Record structure

Flat JSON list, one record per detected atom-atom (or ring-based) contact:

```json
{
  "bgn": {"auth_asym_id", "auth_atom_id", "auth_seq_id", "label_comp_id", "label_comp_type", "pdbx_PDB_ins_code"},
  "end": {same fields as bgn},
  "contact": [list of string labels, non-mutually-exclusive],
  "distance": float (Å),
  "interacting_entities": "INTRA_SELECTION",
  "type": "atom-atom" | "plane-plane" | "atom-plane" | "group-group" | "group-plane"
}
```

## `type` field — observed across the full 533-structure cohort

**100% `"atom-atom"`.** Zero `plane-plane`/`atom-plane`/`group-group`/`group-plane` (ring/π-interaction) records observed across all 533 native VHH structures (confirmed at full-cohort parse time in `stage7_build_tables.py`, not merely the 7/8-structure smoke-test sample). This is treated as a genuine empirical property of this cohort/application (whole native VHH variable domain, no antigen/ligand selection) — the ring-interaction code path is unconditionally invoked (source-confirmed), so this is not a disabled feature. The parser (`stage7_build_tables.py`) handles non-`atom-atom` types generically and would surface them if present in any future re-run on a different cohort.

## `contact` vocabulary (complete possible set, atom-atom records; source: `interactions.py` lines 178-180)

```
clash, covalent, vdw_clash, vdw, proximal, hbond, weak_hbond, xbond, ionic,
metal_complex, aromatic, hydrophobic, carbonyl, polar, weak_polar
```

Multi-label: a single record's `contact` array can hold multiple simultaneous labels (e.g. `["proximal", "weak_polar"]`).

## Classification used in this project's long table (`stage7_interaction_records_long.csv`)

- `is_covalent` = `'covalent' in contact`. **Important real finding:** covers both the trivial CDR3-117/FR4-118 backbone peptide bond AND genuine Cys-Cys disulfide bonds at non-adjacent sequence positions (observed at target positions 50 and 55 in this cohort) — not exclusively a sequence-adjacency artifact.
- `is_proximal_only` = `contact == ['proximal']` exactly (the generic distance-based fallback with no specific chemistry typed).
- `is_specific_interaction` = any label present other than `'proximal'` or `'covalent'`.
- `is_sequential_adjacent` = `|cdr3_imgt_position - target_imgt_position| <= 1` (only ever true at the 117/118 boundary, since all other target positions are far in sequence).
- `include_primary_noncovalent` = `is_specific_interaction AND NOT is_covalent AND target_imgt_position in NONADJACENT_PACKING_TARGET` — this project's primary "real packing interaction" flag.

## Residue identity fields

`auth_asym_id` = single-character surrogate chain (reversible via the per-instance mapping table written alongside each CIF, `data/interim/stage7/cohort_chain_mapping.jsonl`). `auth_seq_id` and `pdbx_PDB_ins_code` are written unchanged from the original structure extraction — no renumbering, confirmed via direct round-trip verification (`stage7_arpeggio_smoke_test.md` §9-10).
