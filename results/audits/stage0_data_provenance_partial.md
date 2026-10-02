> **UPDATE (post Stage-1 kickoff):** FAB+FC has been deliberately removed
> from the primary project scope by user decision. Primary raw cohorts are
> now SD-H, FAB, FV. FAB+FC is retained only as an optional future
> expansion dataset and its absence is no longer a pending-completion item.
> See `results/audits/stage1_metadata_qc.md` for the current status. The
> body of this document is preserved unmodified as the historical Stage-0
> record.

# Stage 0 — Data Provenance Audit (Partial: FAB+FC summary pending)

This audit revalidates the Stage-0 provisional, coordinate-content-based
cohort labels against the official SAbDab2 summary CSVs added to the
project root.

## Official summary CSVs found

| File | Rows | Unique PDB | `type` (exclusive) |
|---|---|---|---|
| `sabdab_search_summary_20260901T161115.csv` | 1597 | 778 | SD-H |
| `sabdab_search_summary_20260901T161202.csv` | 2798 | 1730 | FAB |
| `sabdab_search_summary_20260901T161254.csv` | 536 | 290 | FV |

**No FAB+FC summary CSV is present.** Not assumed, not inferred.

## Matching method

Full extraction manifest (`data/raw_extracted/extraction_manifest.tsv`,
6323 rows, one per original archive member across all 6 `.tgz` archives,
including the duplicate) rebuilt from scratch by re-reading every original
`.tgz` and locating its content on disk — not by trusting the
possibly-renamed extracted filename. Each member's original `(PDB, Hchain,
Lchain)` triple, parsed directly from the *original* tar member path, was
looked up against a combined index of all three summary CSVs keyed on
`(PDB, Hchain, Lchain)`.

## Results

**Match fraction: 100% for every archive, in every direction that matters.**

- SD-H archive (143831): 1597/1597 matched, all `type=SD-H`.
- FV archive (151240): 536/536 matched, all `type=FV`.
- The three previously-"UNKNOWN" Fab-format archives (145541: 373/373,
  150035: 1033/1033, 150654: 1392/1392) **all matched `type=FAB`**, with
  zero unmatched structures. Combined, they account for exactly 2798/2798
  rows of the FAB summary (373+1033+1392=2798), with 0 summary-only rows.
  **This resolves the Stage-0 ambiguity**: these are not a FAB/FAB+FC
  split — they are three non-overlapping, paginated exports of one
  official FAB search. The Stage-0 "UNKNOWN" label is superseded; see
  `archive_summary_crosswalk.tsv` for the authoritative revised mapping.
- The duplicate archive (151141) also matches 1392/1392 FAB rows — the
  *same* 1392 rows already covered by 150654, confirming it is a
  duplicate download, not additional data.

Zero rows were written to `unmatched_downloaded_structures.tsv` — every
downloaded structure found an official metadata match. Zero summary-only
rows were found in any of the three CSVs against this dataset (i.e., in
this case the "broad superset" caveat did not materialize: these three
CSVs happen to align 1:1 with what was downloaded).

## FV archive: confirmed provenance, construct composition not yet filtered

The FV archive is confirmed, via exact official-metadata match (not
coordinate-content inference), to be the genuine output of the SAbDab2 FV
search. Its 536 matched structures break down by the official
`chainsharing_construct` field as:

- 294 — `NA` (conventional non-chain-sharing paired VH/VL, separate
  chains)
- 223 — `SCFV`
- 17 — `DIABODY`
- 1 — `DVD`
- 1 — `CODV`

This confirms and precisely quantifies the Stage-0 coordinate-based
observation ("VH-only + scFv fusion pattern"). Per instruction, these
engineered chain-sharing rows are **not excluded here** — that is a Stage
1/2 scientific-QC decision, not a provenance-audit decision. Provenance
and scientific eligibility are being kept as separate questions.

## FAB archive: official provenance confirmed, not inferred from CH1/CL content

All three FAB-labeled archives are now assigned strictly on the basis of
the official `type=FAB` metadata match, superseding the Stage-0
content-based (CH1/CL presence) reasoning.

## FAB+FC: unresolved, correctly so

No archive is assigned to FAB+FC. No official FAB+FC summary CSV exists in
the project root. RCSB spot-checks and cropped-coordinate inference (used
provisionally in Stage 0) are explicitly **not** used to finalize this
cohort now that official metadata is the standard in use. Status:
`PENDING_OFFICIAL_SUMMARY`.

## Duplicate archive audit

- **Byte-identical .tgz**: none among the 6 archives (all 6 SHA256 values
  differ).
- **Content-identical (decompressed)**: exactly one pair —
  `sabdab_search_structures_20260901T150654.tgz` and
  `sabdab_search_structures_20260901T151141.tgz` — identical sorted tar
  member listing and identical per-member content hashes (differing
  archive-level SHA256 only because gzip embeds a compression timestamp).
  Both matched the same 1392 FAB summary rows, confirming they are the
  same download performed twice. Both original `.tgz` files are retained
  unmodified; only one extracted copy is used downstream (150654).

No other duplicate or content-duplicate archives were found.

## Extraction manifest verification

`data/raw_extracted/extraction_manifest.tsv` rebuilt to cover **every**
member of every archive (6323 rows total = 1597+536+1392+1392+373+1033).
Every original archive member maps to exactly one on-disk extracted file
(0 `MISSING_ON_DISK`). 22 members (18 from the 5 non-duplicate archives +
4 from the duplicate archive's already-covered content) are flagged
`case_collision_flag=True`, meaning they resolve to a
`__casecollision2`-suffixed filename because their exact-case sibling
(`_B_` vs `_b_`, `_H_L_` vs `_h_l_`) already occupied that path on this
case-insensitive filesystem. No crystallographic instance was lost.

## Blocking problems

None found for SD-H, FV, or FAB. FAB+FC provenance is incomplete but this
is expected/documented, not an error, and does not block use of the other
three cohorts.

## Status

**DATA PROVENANCE STATUS: PARTIAL PASS — FAB+FC SUMMARY STILL REQUIRED**
