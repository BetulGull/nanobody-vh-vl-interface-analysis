# Figure 11 — Exemplar Selection and Sequence-Mapping Audit

## Selection rule (deterministic, non-manual)

For each group (VHH, VH) independently:
1. Candidate pool = `is_primary_representative == True` (existing frozen nonredundancy flag, `stage4_vhh_sequence_representatives.parquet` for VHH / `stage3_paired_representatives.parquet` for VH) AND `primary_interface_position_coverage >= 0.90` AND `resolution` present.
2. Standardize (z-score, within that group's candidate pool) 7 common validated structure-level variables (see below).
3. Compute Euclidean distance of each candidate's standardized vector to the pool's own centroid (origin after standardization).
4. Rank ascending; select the 3 nearest-to-centroid structures.

Candidate pool sizes: VHH = 394 (of 533), VH = 743 (of 1429).

## Variables used for the common multivariate selection space

All 7 are present for 100% of both groups (verified by direct null-count check, not assumed):

1. `former_interface_core_PSH`
2. `former_interface_core_PPC`
3. `former_interface_core_PNC`
4. `full_interface_engagement_fraction`
5. `FULL_EMPIRICAL_INTERFACE_28_charge`
6. `FULL_EMPIRICAL_INTERFACE_28_mean_hydropathy`
7. `CDR3_length`

Source: `data/processed/final_figures_model_input.parquet`.

## Variables excluded and why

- `interface_potential_positive_mean` / `interface_potential_negative_mean` (APBS): excluded from the common selection space — not 100% available in both groups (VHH 521/533, VH 1386/1429). Including them would require an additional missing-data branch; the 7 selected variables already give 100%-complete common coverage.
- `shielding_fraction`, `delta_sasa_shield`, `cdr3_sequence_length`, `cdr3_rho`, `cdr3_compactness`, `FR2_CDR3_Dis`, `nonadj_BSA_pair_total`: all VHH-only (0/1429 for VH) — genuinely do not exist for conventional VH, so per the explicit instruction they were never forced into the common space.
- `primary_interface_position_coverage`: used as a candidate-pool quality filter (>=0.90), not as a selection-space phenotype variable, following the coverage-threshold convention already established elsewhere in this project (Figure 5, Figure 6 exemplar selection).

## Selected structures

**VHH** (rank / structure_id / selection_distance):
1. `pdb_00006hjy-G` (0.690)
2. `pdb_00006hjx-F` (0.838)
3. `pdb_00008dqu-A` (0.854)

**Conventional VH**:
1. `pdb_00005ggv-H-L` (0.593)
2. `pdb_00003wfd-H-L` (0.615)
3. `pdb_00007l7d-H-L` (0.681)

3+3 was retained (not reduced to 2+2) — the six-panel grid rendered cleanly with all molecular detail, sequence strips, and IMGT labels legible at the chosen figure size.

## Redundancy handling

All 6 selected structures carry distinct `redundancy_group` values (`sabdab2_H01Q1L0000`, `sabdab2_H01Q0L0000`, `sabdab2_H0437L0000` for VHH; `sabdab2_H0143L00YZ`, `sabdab2_H00NTL00K4`, `sabdab2_H031AL02D8` for VH) — verified directly, not assumed. `pdb_00006hjy-G` and `pdb_00006hjx-F` have very similar phenotype values but are confirmed distinct sequence groups by the existing frozen nonredundancy mapping (not near-duplicate depositions of the same antibody); no additional redundancy logic was invented.

## Framework alignment method

Reused the exact, already-validated common-framework Kabsch alignment machinery (`src/qc/final_manuscript_structures.py::align_to_reference`, `src/qc/stage12_raybould_alignment.py::get_ca_by_imgt/kabsch`, identical to the alignment used in Figures 4A/5/7) — no new alignment logic. The first selected VHH structure (`pdb_00006hjy-G`) was used as the reference frame; the other 5 structures were Kabsch-aligned onto it using common IMGT framework-region CA atoms.

| structure_id | group | RMSD to reference (Å) | N common FR CA |
|---|---|---|---|
| pdb_00006hjy-G | VHH | 0.000 (reference) | 90 |
| pdb_00006hjx-F | VHH | 0.294 | 90 |
| pdb_00008dqu-A | VHH | 0.637 | 90 |
| pdb_00005ggv-H-L | VH | 0.804 | 90 |
| pdb_00003wfd-H-L | VH | 1.403 | 90 |
| pdb_00007l7d-H-L | VH | 1.300 | 90 |

## Important positions mapped

All 14 union-of-top-5 positions from the Figure 10 feature-position matrix (40, 42, 49, 50, 52, 55, 66, 69, 107, 112, 112B, 115, 116, 118) were successfully mapped to an actual amino acid in all 6 selected structures — 0 missing (`Figure11_position_mapping.csv`, `mapped_successfully` all `True`).

## No manual cherry-picking

Every selected structure follows directly and reproducibly from the deterministic centrality rule above, applied to already-frozen structure-level quantities. No structure was chosen, swapped, or excluded based on how its render looked.

## Structures are real experimental entries

All 6 selected structures are genuine SAbDab-derived experimental entries (`INSTANCE` values matching the project's standard `pdb_XXXXXXXX-<chain>` convention, resolved from real deposited coordinate files via `StructureResolver`/`MMCIFStructure`, the same resolution pipeline used throughout this project) — no averaged, synthetic, or median/composite structure was constructed for Figure 11.

## Sequence-visualization revision (color-consistency audit)

Per the follow-up instruction, the original text-list annotation ("IMGT position → feature membership") was removed and replaced with a structure-specific sequence-logo-style strip directly beneath each of the 6 3D panels, using each exemplar's own actual amino-acid sequence (extracted from the same aligned PDB coordinate files used for the 3D renders, chain H only for VH — VL sequence is never included in the sequence strip even though VL remains visible as 3D structural context).

**3D-to-sequence color mapping** (verified identical, not just intended-to-match): the current 3D rendering encodes every important hotspot position with a single gold marker color and the FR2 tetrad (42/49/50/52) with an additional magenta ring — no separate per-feature color scheme was ever rendered in the accepted 3D panels (a distinct per-feature 6-color scheme was tried during development and rejected as visually unreadable given that 5 of the 14 positions carry 4-5 simultaneous feature memberships each). The sequence strips reuse exactly these same two colors: gold fill for every important position, magenta border for tetrad positions, and no highlight for the ±2 context-only neighboring residues. `Figure11_exemplar_sequence_mapping.csv` records, for all 264 (position × structure) rows, both `threeD_color_encoding` and `sequence_color_encoding`; a direct comparison confirms **0 mismatches** between the two columns.

Full per-position feature membership (which of the six Figure-10 features each hotspot belongs to) remains available in `Figure11_position_mapping.csv` and in the Figure 10 source data / caption — it is documentation, not a repeated on-figure color code, since the accepted 3D rendering itself does not encode individual features by separate colors.

## Local sequence window construction

For each of the 14 hotspot positions, the ±2 immediate sequence neighbors (true chain order, read directly from each structure's own resolved residue list — not a naive IMGT-integer ±2, which would be wrong across CDR3 insertion-code boundaries) were included as neutral context. Overlapping/adjacent windows were merged; non-adjacent windows are rendered as visually separate blocks with a dashed divider, grouped as FR2 / FR3 / CDR3 / FR4, so that widely separated true-sequence gaps (e.g. position 57 to position 64) are never shown as if contiguous.

## No new scientific calculations

No interface derivation, FreeSASA, shielding, antigen-contact, APBS, patch, or structure-selection analysis was rerun. The only new computation in this entire figure is the deterministic multivariate-centrality exemplar ranking (a visualization-selection procedure, not a biological metric) and a direct read of existing amino-acid identities from already-resolved, already-aligned coordinate files.
