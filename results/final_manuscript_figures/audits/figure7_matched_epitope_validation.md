# Figure 7 — Matched-Epitope Validation Gate

Candidate pair identified by an automated antigen-name scan across the frozen Stage-9 contact table (`src/qc/final_freeze_3d_exemplars.py`): **VHH = `pdb_00008k4q-B`**, **VH = `pdb_00006wgl-A-B`**, reported antigen "Interleukin-4 receptor subunit alpha". Raw 4Å antigen-contact Jaccard (own residue numbering, unremapped) = 0.92 (11/12 shared). This candidate is treated as unproven until the steps below pass.

## 1-2. Antigen chain identification and sequence extraction

Identified directly from the frozen per-instance mmCIF files via `gemmi` polymer entity extraction (not assumed from metadata):
- `pdb_00008k4q-B`: antigen = **chain A**, 198 polymer residues (author numbering 3-200), sequence starts `VLQEPTCVSDY...`
- `pdb_00006wgl-A-B`: antigen = **chain C**, 178 polymer residues (author numbering 1-196), sequence starts `MKVLQEPTCVSDY...`

Both confirmed as the sole antigen chain per instance against the frozen Stage-9 `antibody_antigen_contacts_long` table (`antigen_chain` == 'A' for 8k4q, == 'C' for 6wgl, no other antigen chain present in either).

## 3. Sequence alignment

Global pairwise alignment (Biopython `Align.PairwiseAligner`, match=+2/mismatch=-1/gap-open=-10/gap-extend=-0.5) of the two extracted polymer sequences.

## 4. Alignment statistics

| metric | value |
|---|---|
| VHH-complex antigen chain length | 198 residues |
| VH-complex antigen chain length | 178 residues |
| Aligned (matched, ungapped) length | 176 residues |
| Sequence identity over aligned region | **99.43%** (175/176) |
| Coverage vs VHH-complex antigen | 88.9% |
| Coverage vs VH-complex antigen | 98.8% |
| Mismatches in aligned region | 1 (8k4q residue 182 CYS vs 6wgl residue ~164 LEU — outside the contact region, likely a construct/allele difference) |
| Indels | 6wgl chain C carries an N-terminal `MK` extension (expression-tag/signal remnant) absent from 8k4q chain A's modeled range; three short internal gaps (4, 4, 9 residues) in disordered/unresolved loop stretches — handled natively by the aligner's affine gap penalty, not by any manual truncation |

This is essentially the same IL-4Rα ectodomain construct solved in two independent depositions — not two different antigens sharing a name by coincidence.

## 5. Residue mapping

Full 1:1 residue mapping (author-numbered) written to [`figure7_antigen_residue_alignment_mapping.csv`](../tables/figure7_antigen_residue_alignment_mapping.csv) (176 rows: pdb1_resnum/pdb1_aa ↔ pdb2_resnum/pdb2_aa, `match` flag).

## 6-7. Remapped contact sets

The VH-complex (6wgl) contact set was remapped into the VHH-complex (8k4q) author-numbering coordinate system via the table above, then intersection/union/Jaccard recomputed.

| | raw numbering | remapped |
|---|---|---|
| VHH contacts (4Å) | {41,42,43,44,45,47,66,67,69,71,72} (n=11) | unchanged (reference frame) |
| VH contacts (4Å) | {38,41,42,43,44,45,47,66,67,69,71,72} (n=12) | identical set — every VH contact residue number mapped 1:1 to the same number in 8k4q's frame (n=0 unmapped) |
| Intersection | — | 11 |
| Union | — | 12 |
| Jaccard | 0.9167 (raw) | **0.9167 (remapped)** |

## 8. Numbering-artifact check

**No artifact.** The remapped Jaccard is numerically identical to the raw-numbering Jaccard (0.9167 both ways) because the two depositions already share the same author-numbering convention for this construct over the entire contact region — confirmed independently via the alignment (every one of the 12 relevant residue numbers falls in a 1:1-matched, identical-residue-name aligned block; residue names checked directly: e.g. both structures have PHE41, LEU42, LEU43, SER44, GLU45, HIS47, ASP66, ASP67, VAL69, ALA71, ASP72). The apparent overlap is real, not a coincidence of independent numbering schemes colliding.

## 9. Spatial co-localization (structural superposition)

Independent structural confirmation, not just a sequence-numbering argument: Kabsch superposition of the two antigen chains on all 176 aligned Cα pairs gives a whole-chain RMSD of **2.34 Å** (176 CA pairs) — reasonable for two independently solved crystal structures of the same ectodomain, not evidence of misassignment. Restricting to the 11 shared contact residues after applying that same whole-chain superposition gives an RMSD of **2.11 Å** (11 CA pairs) — essentially the same magnitude as the whole-chain RMSD, confirming the shared contact residues genuinely co-localize in 3D space after alignment, not merely by coincident residue numbers.

## Verdict

**Status: `VALID_MATCHED_EPITOPE`**

All required checks pass: same antigen construct (99.4% identity), contact-residue overlap confirmed both in the remapped numbering (Jaccard=0.917) and independently in 3D space after structural superposition (epitope-subset RMSD ≈ whole-chain RMSD). Figure 7 Panels G/H may therefore describe this pair as "matched antigen / overlapping epitope" — this is a genuine, validated finding for this cohort, not assumed by analogy to Gordon et al. 2023 Figure 14.

Full machine-readable record: [`figure7_matched_epitope_validation.csv`](../tables/figure7_matched_epitope_validation.csv).
