# Stage 7 — Covalent-Adjacency Target-Set Policy

**Derivation:** from actual IMGT sequence order (not assumption), cross-referenced against the frozen Stage-3 empirical interface (`stage3_empirical_vh_vl_interface_primary.csv`, 28 primary positions) and the Stage-6 CDR3 definition (IMGT 105-117, `results/audits/stage6_gordon_method_definition.md` / shielding audit).

## 1. Frozen inputs

- Empirical interface (28 positions, verified count from file, not hard-coded):
  `40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 111A, 112, 112A, 112B, 113, 114, 115, 116, 118, 119, 120`
- CDR3 (IMGT 105-117, including any insertion codes in that numeric range): the interface positions falling inside this range are
  `107, 108, 109, 110, 111, 111A, 112, 112A, 112B, 113, 114, 115, 116` — **13 positions**, matching Stage 6's frozen "13 of 28 fall within CDR3" figure exactly.
- **FULL_SHIELDABLE_TARGET** (28 − 13 = 15 positions, all in FR2/FR3/FR4, reused unchanged from Stage 6):
  `40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120`

## 2. Covalent-adjacency determination

IMGT numbering places CDR3's defined final position at 117 (insertions in this dataset occur only at 111/112, never after 117), and FR4 begins immediately at 118 with no intervening real residue — i.e., 118 is, by construction of the IMGT numbering scheme itself (not an assumption about any specific structure), always the immediate sequence successor of the last CDR3 residue. **118 is therefore directly covalently bonded to CDR3's C-terminal residue** and is excluded from the primary packing target set.

On the N-terminal side, CDR3 begins at IMGT 105; its immediate sequence predecessor is IMGT 104, which is **not** a member of the empirical interface (only 103 is, and 103 is separated from CDR3 by 104 — i.e., 103 is not covalently adjacent to CDR3, it is one residue removed). No N-terminal covalent-adjacency exclusion applies.

This was verified directly from the interface CSV's own `imgt_position`/`insertion_code` columns and the CDR3 boundary definition already frozen in Stage 6 — not inferred from memory of IMGT conventions alone.

## 3. Target sets (frozen for Stage 7)

- **NONADJACENT_PACKING_TARGET** (primary analysis) = FULL_SHIELDABLE_TARGET − {118} =
  `40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120` — **14 positions**
- **FULL_SHIELDABLE_TARGET** (required sensitivity analysis, unchanged from Stage 6) =
  `40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120` — **15 positions**

Both sets are used throughout Stage 7: NONADJACENT_PACKING_TARGET as the primary BSA/interaction analysis (avoiding conflating genuine noncovalent packing with the trivial geometric consequence of a peptide bond at 118), FULL_SHIELDABLE_TARGET reported in parallel as a sensitivity check.
