# Figure 4A/4B — Input Lineage Audit

Read-only investigation completed before any new computation. All findings verified directly against frozen files (not from memory).

**Naming note**: this project's `results/final_manuscript_figures/` already contains an earlier, separate "Figure 4" ("FR2-CDR3 geometry", VHH-only quantile exemplars) and content in "Figure 3"/"Figure 6" that touches CDR3 length. Per this task's explicit scope restriction, that prior figure is **out of scope, untouched, not referenced further**, and not registered under the new `Figure4A`/`Figure4B` manifest entries. The new outputs use distinct filenames (`Figure4A.*`, `Figure4B.*`) and never overwrite `Figure4.*`.

## 1. Position-set definitions (verified exact)

| set | size | positions | source |
|---|---|---|---|
| `FULL_EMPIRICAL_INTERFACE` (code constant) | 25 unique integers | `{40,42,44,49,50,52,55,66,67,68,69,103,107,108,109,110,111,112,113,114,115,116,118,119,120}` | `src/qc/final_regional_physicochemistry.py`, `stage8_apbs_full_cohort.py`, `stage8_patch_full_cohort.py`, `stage8_topology_full_cohort.py`, `stage12_raybould_write_structures.py` — identical literal across all |
| Stage-3 "28-position" table | 28 rows | same 25 integers + 3 extra insertion-coded CDR3 sub-positions (111A, 112A, 112B) | `results/tables/stage3_primary_interface_positions.csv` — a row-count label, not a different position set (see the Figure-1 25-vs-28 audit for the full reconciliation) |
| `SHIELDABLE_NON_CDR3_INTERFACE` | **15** | `{40,42,44,49,50,52,55,66,67,68,69,103,118,119,120}` | `src/qc/stage6_cdr3_removal.py::get_shieldable_interface_positions()` = `FULL_EMPIRICAL_INTERFACE` minus CDR3 range (105-117); independently confirmed identical to `stage7_bsa.py`'s `FULL_SHIELDABLE_TARGET` |
| `NONADJACENT_PACKING_TARGET` | **14** | `{40,42,44,49,50,52,55,66,67,68,69,103,119,120}` | `SHIELDABLE_NON_CDR3_INTERFACE` minus position 118 (covalently adjacent to CDR3's C-terminal end; policy in `results/audits/stage7_covalent_adjacent_policy.md`) |

Figure 4A's Panel B/C/D/E general geometry map uses **`SHIELDABLE_NON_CDR3_INTERFACE` (15 positions)**, not the 14-position packing target, per the task's explicit instruction (Section 3). The 14-position set is reserved for Stage-7-style packing-burden calculations and is not used to redefine Figure 4A's target.

## 2. Existing CDR3-geometry data — coverage gap

- `results/tables/stage7_interaction_records_long.csv` (89,725 rows) is **VHH-only** (verified: `INSTANCE.nunique()==533`, zero VH rows) and is a **discrete contact list, not a continuous distance matrix** — every row's `distance` is an Arpeggio-detected contact within Arpeggio's own implicit ~5.0 Å detection radius (confirmed empirically: `distance.max() == 5.00` exactly across all rows). Residue pairs farther apart simply have no row. This table cannot supply a VHH-vs-VH comparison and cannot supply true minimum-distance geometry beyond ~5 Å — a **fresh computation is required for both groups**.
- `data/processed/final_figures_model_input.parquet`: `FR2_CDR3_Dis`, `cdr3_rho`, `cdr3_compactness`, `nonadj_BSA_pair_total` are **VHH N=533, VH N=0** (100% NaN for VH). Only `CDR3_length` and `local_thr4.0_largest_patch_residue_count` are populated for both groups. Confirms: any VHH-vs-VH comparison of CDR3 geometry, packing distance, or FR2_CDR3_Dis requires new VH-side computation — none of it is frozen for VH.

## 3. FR2_CDR3_Dis — frozen definition, reused exactly

`src/qc/stage6_geometry.py::compute_zhang_fr2_cdr3_dis()`: centroid of all non-H, altloc-resolved atoms at IMGT 42, minimum Euclidean distance to the centroid of each internal CDR3 residue (real CDR3 residues with the first and last two excluded, i.e. `cdr3_residues[2:-2]` in true sequence order). Matches Zhang et al. 2025 Methods §2.2 verbatim (independently re-extracted from the source PDF in this audit, see `figure4B_zhang_comparative_metric_audit.md`). VHH-only in the frozen table; the same function is reused unmodified for the new VH-side calculation in this task (Section 5 below), applied to native VH coordinates from the paired VH-VL structure with VL simply never loaded (not "isolated," not remodeled).

## 4. DSSP

No existing DSSP output anywhere in the repo for either group. `mkdssp` (v4.6.1) is installed and functional in the `vhh3d` environment (confirmed via `stage0_environment_setup.md` and a live version check) — a fresh, lightweight DSSP run on native experimental coordinates (no minimization/relaxation) is required for Figure 4B Panel E's secondary-structure metrics.

## 5. Structure-extraction infrastructure — reused, not rebuilt

`src/qc/final_manuscript_structures.py` (built in the prior final-manuscript-figures session) already provides `get_variable_domain_atoms(instance, chain_which)` (IMGT-numbered, altloc-resolved, no-H atoms for either the H or L chain of any instance in the combined VHH+VH representative lookup), `align_to_reference`, and `write_pdb`. Calling `get_variable_domain_atoms(instance, 'H')` on a VH instance naturally extracts only heavy-chain atoms — VL is a physically separate chain and is never loaded, satisfying "VL omitted from the mathematical surface/geometry calculation; VH coordinates retained unchanged from the native paired Fv" without any special-casing. This module is reused directly for all new atom-level computation in Figure 4A and 4B.

## 6. Contact-distance threshold convention

The project's own existing, declared, "Gordon-verified primary" heavy-atom contact threshold is **4.5 Å** (`src/qc/contact_search.py`, `stage9_antigen_contacts.py`), with 4.0/5.0 Å as declared sensitivity thresholds (`stage9_contact_threshold_sensitivity.csv`, both groups, all three thresholds already validated there). **Figure 4A adopts this exact convention: primary=4.5 Å, sensitivity={4.0, 5.0} Å** — an existing, validated project convention, not a new arbitrary choice, and not selected by which result looks strongest (declared here, before any contact-probability values are computed).

Zhang et al. 2025's own interaction threshold is 4 Å (their Methods, re-extracted independently — see the Zhang metric audit). Figure 4B's CDR3-CDR1/FR2/CDR2 interaction panels use Zhang's own 4 Å threshold to stay faithful to the reproduced analytical logic, explicitly distinguished in the caption from Figure 4A's project-native 4.5 Å convention.

## 7. Zhang et al. 2025 exact-metric recovery

See the dedicated companion audit: [`figure4B_zhang_comparative_metric_audit.md`](figure4B_zhang_comparative_metric_audit.md).

## 8. Summary of what is newly computed vs reused

| quantity | status |
|---|---|
| 15/14/25-position sets | reused exactly, no redefinition |
| FR2_CDR3_Dis (VHH) | reused frozen value |
| FR2_CDR3_Dis (VH) | **new**, same frozen function, new coordinates |
| CDR3-to-15-position minimum heavy-atom distance (both groups) | **new** — no equivalent frozen table exists |
| CDR3-to-CDR1/FR2/CDR2 minimum heavy-atom distance (both groups) | **new** |
| CDR3 secondary structure (DSSP G/T/E/S/C) | **new**, both groups |
| Domain SASA, SASA×hydrophobicity ("Zhang surface hydrophobicity") | **new**, both groups, VL mathematically omitted for VH |
| Zhang LPP/LNP (PEP-Patch-derived) | **NOT reproduced** — PEP-Patch is not part of this project's toolchain; omitted from the exact Panel-E metric set per the task's explicit instruction not to silently substitute |
| Zhang "surface charge" (PEP-Patch-derived patch feature) | **NOT reproduced** for the same reason; the project's own residue-level net-charge proxy is a different quantity and is not relabeled as Zhang's metric |
| Solubility | **not computed**, intentionally, per explicit instruction |
