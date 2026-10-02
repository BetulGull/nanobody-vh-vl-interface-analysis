# Figure 6 — Authoritative BSA Definition

Source: `src/qc/stage7_bsa.py` (frozen, unmodified, reused verbatim for both groups in this revision).

## Residue sets

- **CDR3** (side A): IMGT 105-117 inclusive, ALL real (observed) residues with any insertion codes present — `CDR3_START, CDR3_END = 105, 117`. Not restricted to "internal" residues (that exclusion is specific to the FR2_CDR3_Dis metric elsewhere in this project, not to BSA).
- **Packing target** (side B): the frozen 14-position `NONADJACENT_PACKING_TARGET` = `{40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}`. This is the 15-position `SHIELDABLE_NON_CDR3_INTERFACE` minus position 118 (excluded because it is covalently adjacent to CDR3's C-terminal end — `results/audits/stage7_covalent_adjacent_policy.md`). **This exact 14-position set is preserved for the VH calculation — never substituted with the 15-position set.**

## SASA parameters

- Engine: FreeSASA (freesasa-python 2.2.1), `src/qc/stage3_sasa.py::calc_sasa`, the same validated engine used throughout Stages 3/5/6.
- Algorithm: Lee & Richards (explicit).
- Probe radius: 1.4 Å (explicit).
- Classifier: FreeSASA default (NACCESS-derived van der Waals radii).
- Atoms: all non-hydrogen atoms (backbone + side chain both included), altloc-resolved (highest occupancy, alphabetical tie-break), from `stage3_variable_domain.extract_variable_domain_atoms` (IMGT 1-128 only).

## Exact formula

```
SASA_A_alone   = total SASA of CDR3 atoms, computed in isolation (CDR3 only, no target atoms present)
SASA_B_alone   = total SASA of packing-target atoms, computed in isolation (target only, no CDR3 atoms present)
SASA_AB_total  = total SASA of CDR3 + packing-target atoms computed TOGETHER, in native coordinates (no other atoms — not the whole domain, not VL, not antigen)

BSA_pair_total = SASA_A_alone + SASA_B_alone - SASA_AB_total
```

This is the standard buried-surface-area convention (sum of the two isolated areas minus the area of the combined system) — **not divided by 2**. A per-side additive decomposition is also computed and cross-validated:

```
BSA_CDR3_side   = SASA_A_alone - (sum of per-residue SASA over CDR3 residues, evaluated within the AB complex)
BSA_target_side = SASA_B_alone - (sum of per-residue SASA over target residues, evaluated within the AB complex)
```

with the required consistency check `BSA_CDR3_side + BSA_target_side == BSA_pair_total` (exact by construction; verified `< 1e-6` tolerance for every structure in both groups in this revision).

## Missing-residue handling

If either the CDR3 atom set or the target atom set is empty for a given structure (e.g. a target position entirely absent/unresolved), the structure is marked `status='TARGET_NOT_EVALUABLE'` and excluded from the BSA computation — never zero-filled.

## Coordinate context

- **VHH**: single heavy-chain coordinates (no VL exists in a VHH structure).
- **Conventional VH (this revision)**: heavy-chain-only atoms extracted from the native, paired VH-VL experimental structure via `get_variable_domain_atoms(instance, 'H')` — VL is a physically separate chain and is simply never loaded into the atom list passed to `compute_pairwise_bsa`. This matches the function's own existing "single heavy-chain coordinate context" input exactly, with no VH remodeling, relaxation, or energy minimization, and no VL geometric deletion from the native structure file (VL coordinates remain on disk, untouched; they are only absent from this specific in-memory atom list). Documented per the task's explicit instruction: this is "VL omitted from the calculation," never described as an experimentally isolated VH domain.
