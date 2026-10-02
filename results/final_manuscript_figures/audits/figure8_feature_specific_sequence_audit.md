# Figure 8 — Feature-Specific Local Sequence Audit

## Authoritative source for important positions (per feature)

All six top-5 important-position lists are read directly from the already-validated, already-audited `results/final_manuscript_figures/source_data/Figure10/Figure10_feature_position_matrix.csv` (`important==True` rows) — the same predeclared, non-manual `TOPN=5` ranking rule documented in `figure10_feature_position_audit.md`. No ranking was recomputed or re-derived for Figure 8; this figure only re-visualizes the already-frozen selection at the sequence level.

| Feature | Important positions (verified from Figure10 source data) |
|---|---|
| Sequence divergence | 40, 42, 49, 50, 52 |
| Hydropathy shift | 42, 49, 50, 52, 115 |
| Charge shift | 49, 50, 69, 107, 116 |
| Native exposure difference | 49, 50, 52, 112, 112B |
| CDR3 shielding (VHH only) | 42, 52, 55, 66, 118 |
| Antigen engagement | 42, 49, 50, 52, 118 |

These exactly match the values previously reported in this workflow (with the trivial exception that Figure 10's own printed values used `abs()`-ranking for signed features, which is what produced this same list — reused verbatim, not recomputed).

## Local sequence-window construction

For each important position, up to 2 true-sequence neighbors on each side were included as context, using one of two rules (documented explicitly to avoid inventing cohort-length-dependent insertion-code neighbors):

- **Plain-integer positions** (all except the two CDR3 insertion-anchor positions 112 and 112B): true integer ±1/±2 neighbors (e.g. position 40 → context 38, 39, 41, 42).
- **112 and 112B** (the two positions that themselves carry/adjoin cohort-variable insertion codes): context drawn from the already-validated 28-entry Stage-3 canonical interface order (`111, 111A, 112, 112A, 112B, 113`) rather than a true single-structure insertion-code neighbor, which would differ by CDR3 length across the cohort and cannot be defined at the cohort-logo level. This is a deliberate, documented simplification, not an invented ranking.

Overlapping/adjacent windows were merged into contiguous blocks; non-adjacent windows are rendered as visually separate local blocks with a clear gap, so that distant true-sequence positions are never drawn as if adjacent (e.g. Charge shift's FR2 block [47-52], FR3 block [67-71], and two separate CDR3 blocks [105-109] and [114-118] are four visually distinct blocks, never one continuous strip).

## Local blocks shown per feature

| Feature | Blocks (region: positions) |
|---|---|
| A. Sequence divergence | FR2: 38-44 \| FR2: 47-54 |
| B. Hydropathy shift | FR2: 40-44 \| FR2: 47-54 \| CDR3: 113-117 |
| C. Charge shift | FR2: 47-52 \| FR3: 67-71 \| CDR3: 105-109 \| CDR3: 114-118 |
| D. Native exposure | FR2: 47-54 \| CDR3: 111,111A,112,112A,112B,113,114 |
| E. CDR3 shielding (VHH only) | FR2: 40-44 \| FR2: 50-57 \| FR3: 64-68 \| FR4: 116-120 |
| F. Antigen engagement | FR2: 40-44 \| FR2: 47-54 \| FR4: 116-120 |

## Sequence-logo data source

Position-level VHH and conventional-VH amino-acid frequencies for every displayed position (44 unique positions across all six panels) were read directly from the existing frozen `data/processed/stage4_position_aa_frequencies.csv` (the same table underlying the Figure 1/2/9 sequence logos), X-dropped and renormalized per (group, position) — the identical convention already established and audited for those earlier logos. No new sequence curation, alignment, or frequency computation was performed.

## Region labels

FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4 boundaries use the exact IMGT ranges already established and used elsewhere in this project (`scripts/figure4A_render_structural_context.py` / `figure4B_render_panelA.py`: CDR1=27-38, FR2=39-55, CDR2=56-65, FR3=66-104, CDR3=105-117, FR4=118-128).

## Group labels

Logos are labeled "VHH" and "VH" (never "VH-VL") — both logos are heavy-variable-domain sequence summaries only; light-chain sequence is not part of either logo, consistent with the empirical former-interface being a VH-side positional definition.

## Tetrad annotation discipline

The canonical FR2 tetrad (42/49/50/52) is labeled in bold gold on the IMGT position axis in every panel where it appears (a constant structural fact, independent of which feature panel is shown), but only receives the magenta highlight *outline* when it is also an important hotspot for that specific feature (e.g., position 52 appears in panel C, Charge shift, as plain context — gold label, no magenta ring — because 52 is not one of the 5 charge-shift hotspots there).

## No new biological metric invented

This figure performs only: (1) reading already-validated important-position lists from Figure 10's frozen source data, (2) reading already-validated cohort amino-acid frequencies from Stage 4, (3) lightweight windowing/merging/plotting for visualization. No interface derivation, sequence curation, FreeSASA, Arpeggio, APBS, PEP-Patch/TPatch, or antigen-contact calculation was rerun. No solubility, aggregation propensity, stability, or affinity/KD/kon/koff value was calculated anywhere in this figure.
