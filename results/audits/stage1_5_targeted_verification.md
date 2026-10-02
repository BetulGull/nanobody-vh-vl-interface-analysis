# Stage 1.5 — Targeted Provenance / Chain Verification (read-only)

No cohort membership, `sdh_classification`, or `primary_vhh_candidate_metadata`
value was changed. No SASA/geometry/interaction/electrostatics/aggregation/
statistics were calculated. All evidence below comes from actual mmCIF
entity/polymer/coordinate data read via `gemmi`, plus the existing Stage-1
metadata table.

## 1–2. The five conflicting SABDAB_IDs: actual sequence comparison

For every instance under each of the five IDs, the SabDab2 Hchain's
coordinate-observed sequence was extracted with its residue numbering
(already IMGT-numbered by SAbDab2 per project convention), and compared
across instances **only at residue-number positions present in every
compared structure** (so unresolved/missing residues never masquerade as
mismatches).

| SABDAB_ID | n instances | n unique full sequences | n shared IMGT positions compared | mismatches at shared positions | identical at shared positions |
|---|---|---|---|---|---|
| H00QEL0000 | 10 | 4 | 110 | 0 | **True** |
| H021EL0000 | 12 | 2 | 115 | 0 | **True** |
| H02DAL0000 | 2 | 2 | 128 | 0 | **True** |
| H02Y5L0000 | 11 | 3 | 112 | 0 | **True** |
| H06BKL0000 | 2 | 1 | 125 | 0 | **True** |

**Result: zero amino-acid mismatches at any shared IMGT position, in all
five IDs**, despite `n_unique_entity_full_sequences` > 1 in four of them.
The "unique full sequence" count differs only because the full
`entity_poly` sequence includes construct-specific flanking material
(His-tags, short linker residues, minor C-terminal truncation/extension)
that varies by expression construct — the actual IMGT-numbered variable
domain itself is identical everywhere it can be compared. Full detail:
`results/audits/stage1_5_sequence_identity_check.csv`.

This is direct structural evidence, not an assumption, and it resolves in
favor of outcome (1) from the four possible outcomes: **identical
sequence, conflicting provenance/species annotation** — not a genuine
sequence difference (outcome 3) and not a mapping/parsing error (outcome
4, ruled out by the successful extraction and consistent numbering ranges
across every instance).

## 3. Provenance-conflict resolution (no rule/data changes applied)

`results/audits/stage1_5_sabdab_id_provenance_conflicts.csv`:

| SABDAB_ID | species text labels | heavy_taxid values | self-contradictory rows | status |
|---|---|---|---|---|
| H00QEL0000 | homo sapiens; lama glama | 9606; 9844 | 0 | MIXED_PROVENANCE_SAME_SEQUENCE |
| H021EL0000 | homo sapiens; lama glama | 9606; 9844 | 0 | MIXED_PROVENANCE_SAME_SEQUENCE |
| H02DAL0000 | homo sapiens; lama glama | 9606; 9844 | 0 | MIXED_PROVENANCE_SAME_SEQUENCE |
| H02Y5L0000 | homo sapiens; lama glama | 9606; 9844 | 0 | MIXED_PROVENANCE_SAME_SEQUENCE |
| H06BKL0000 | vicugna pacos (both rows) | 30538; 9606 | **2** | **POSSIBLE_DATABASE_INCONSISTENCY** |

Four of the five IDs (H00QE, H021E, H02DA, H02Y5) show the same identical
variable-domain sequence carrying different `heavy_taxid` annotations
(camelid vs. `9606`=Homo sapiens) across different PDB entries — a
cross-entry annotation inconsistency, not evidence the underlying
molecule changed identity.

**H06BKL0000 is a distinct and more specific finding**: in this case
`heavy_species` literally says `"vicugna pacos"` (text) while
`heavy_taxid` says `9606` (number) **within the same row**
(`pdb_00009q00`), for what the `compound` field confirms is the *same*
construct ("Genetically encoded fluorescent nanobody H11") as its sibling
entry `pdb_00009pzz`, whose row correctly has `heavy_taxid=30538`. This is
a self-contradictory record in the source SAbDab2 metadata itself — not
merely a cross-format labeling difference — and is flagged
`POSSIBLE_DATABASE_INCONSISTENCY` accordingly. Note the Stage-1 classifier,
by design, gives `heavy_taxid` top priority in the evidence hierarchy, so
this specific row (`pdb_00009q00`) is currently (incorrectly, on the
weight of this new evidence) routed to `human_sdVH_candidate` even though
`heavy_species` text and `compound` both independently agree it is
Vicugna pacos nanobody H11.

## 4–5. H02Y5L0000: is the MHC single-chain trimer the antibody chain, or the antigen?

Checked at the mmCIF entity/chain level for **all 11 instances** (not a
sample), using actual sequence content — these files carry no
`_entity.pdbx_description` field at all, so identity was determined from
the sequences themselves (a distinctive 40-residue substring of human
beta-2-microglobulin, and the conserved HLA class-I mature N-terminus
`SHSMRYF`).

**Worked example — `pdb_00007sr0`** (representative of the cluster):

| Chain | Length | Content |
|---|---|---|
| A | 424 aa | HPV-16 E7 peptide + GS-linker + beta-2-microglobulin + GS-linker + HLA-A\*02:01 heavy chain (ECD) — a single-chain pMHC-I trimer, C-terminal His-tag |
| **B (SabDab2 Hchain)** | **115 aa** | `EVKLVESGGGLVQPGGSLRLSCAASGSIFSINTMGWYRQTPGKQRDLVADISSGGSTKYGDSVKGRFTISRDNTKNTVYLQMNSLKPEDTAVYYCYGLSYSNDDYWGQGTQVTVSS` |

Chain B is a single, unfused, ~115-residue Ig variable domain: canonical
VH-family FR1 (`EVKLVESGGGLVQPGGSLRLSCAAS...`), FR2 opening
`WYRQTPGKQRDLVA` (qualitatively the classic camelid FR2 substitution
pattern relative to conventional VH's `WVRQ...`, reported here only as
supporting information, not as the basis for the identification — per
instruction 6 IMGT 42/49/50/52 residue identity was not used to decide
this), a single CDR3 (`YGLSYSNDDY`), and FR4 ending `WGQGTQVTVSS`. There
is no internal GS-linker, no second Ig domain, and no MHC/B2M sequence
anywhere within chain B.

**Verified across all 11 instances** (`results/audits/stage1_5_h02y5_chain_verification.csv`):

| PDB | Hchain | Hchain length | Hchain is MHC construct | Hchain looks like VHH | MHC/B2M found in |
|---|---|---|---|---|---|
| pdb_00007sqp | B | 116 | False | True | chain A (427 aa, B2M signature) |
| pdb_00007sr0 | B | 115 | False | True | chain A (424 aa, B2M signature) |
| pdb_00007sr3 | B | 112 | False | True | chain A (427 aa, B2M signature) |
| pdb_00007sr4 | B | 112 | False | True | chain A (429 aa, B2M signature) |
| pdb_00007sr5 | B | 113 | False | True | chain A (424 aa, B2M signature) |
| pdb_00007srk | B | 115 | False | True | chain A (423 aa, B2M signature) |
| pdb_00007ssh | B | 115 | False | True | chain A (429 aa, B2M signature) |
| pdb_00007st3 | B | 115 | False | True | chain A (429 aa, B2M signature) |
| pdb_00007stg | B | 115 | False | True | chain A (423 aa, B2M signature) |
| pdb_00009min | E | 116 | False | True | chain B (101 aa, B2M signature) |
| pdb_00009o5s | F | 117 | False | True | chain D (99 aa, B2M signature) |

**Answer: the MHC single-chain trimer (or, in `pdb_00009min`/`pdb_00009o5s`,
the separately-chained B2M/HLA components alongside a de-novo-designed
minibinder antigen) is a separate antigen/complex component. It is never
the SabDab2 Hchain.** The SabDab2 Hchain is, in every one of the 11
instances, the same short genuine VHH-format sequence. `Lchain` is `NA`
in all 11, as expected for SD-H.

## 6. IMGT hallmark note (supporting only, not decisive)

Per instruction, IMGT 42/49/50/52 identity was not used to decide
anything above — the decision rests on SabDab2 chain assignment, actual
sequence content, entity separation, and provenance metadata. The FR2
pattern noted above is offered only as a qualitative consistency check.

## 7. Reassessment of the "58 contaminant rows" claim

**FALSE_POSITIVE_FROM_COMPOUND_TEXT.** The Stage-1 audit's compound-text
keyword sweep ("single chain trimer", "HLA", etc.) correctly found the
vocabulary, but incorrectly attributed it to the antibody chain rather
than the antigen. Direct entity/chain verification across all 11 PDBs
(58 instance rows) shows the SabDab2 Hchain is a genuine, structurally
and sequence-distinct VHH domain in every case; the MHC/HLA/B2M material
belongs to a separate antigen entity. **None of the 58 rows are
demonstrated non-antibody contamination.** This is exactly the failure
mode the user's correction 2 warned about.

## 8. Conservative policy for mixed-provenance IDs (recommended, not applied)

For the four IDs confirmed `MIXED_PROVENANCE_SAME_SEQUENCE`
(H00QEL0000, H021EL0000, H02DAL0000, H02Y5L0000) and the one
`POSSIBLE_DATABASE_INCONSISTENCY` (H06BKL0000):

- Set `sequence_provenance_ambiguity = true` for every instance under
  these 5 SABDAB_IDs (108 rows) — **not applied in this audit**, listed
  only as a recommendation.
- For later physicochemical deduplication (one representative per
  SABDAB_ID), do **not** select a representative from these 5 IDs without
  first resolving which `heavy_taxid` annotation is correct — treating
  the differently-labeled instances as independent sequences would be
  wrong, since the actual variable-domain sequence is identical.
- Recommend, as a separate follow-up, escalating the `pdb_00009q00`
  self-contradictory row (H06BKL0000) as a likely genuine SAbDab2 database
  error (heavy_species text vs. heavy_taxid number disagree within one
  row) for correction or explicit override, since here the internal
  evidence (species text + compound text) both agree and only one field
  (heavy_taxid) is the outlier.
- No exclusion or reclassification is applied at this stage.
