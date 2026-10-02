# Figure 10 (rebuilt) — Sequence-to-Property Audit

This is a full visual rebuild of Figure 10, replacing the earlier tile-matrix-only version. The underlying important-position rankings, sequence frequencies, and per-position feature values are all **reused unmodified** from already-validated, already-audited sources (the prior Figure 10 tile-matrix build and the Figure 8 local-sequence-window build) — nothing was recomputed.

## Every highlighted position is validated

All "important" positions in every module are read directly from `Figure10_feature_position_matrix.csv::important==True` (the frozen, predeclared top-5-by-value ranking already audited in the original `figure10_feature_position_audit.md`). No position was added, removed, or re-ranked for this rebuild.

| Feature | Important positions |
|---|---|
| A. Sequence divergence | 40, 42, 49, 50, 52 |
| B. Hydropathy shift | 42, 49, 50, 52, 115 |
| C. Charge shift | 49, 50, 69, 107, 116 |
| D. Native exposure | 49, 50, 52, 112, 112B |
| E. CDR3 shielding (VHH only) | 42, 52, 55, 66, 118 |
| F. Antigen engagement | 42, 49, 50, 52, 118 |

## Local sequence windows are correct

Each local window = the important position plus up to 2 true-sequence neighbors each side (plain integer neighbors for all positions except the two CDR3 insertion-anchors 112/112B, which use the already-validated 28-entry canonical order instead of an undefined, cohort-length-dependent insertion-code neighbor — identical rule already documented and reused from the Figure 8 audit). Overlapping windows are merged into one block; non-adjacent windows are rendered as visually separate local blocks with a labeled gap (e.g. Charge shift: FR2[47-52] | FR3[67-71] | CDR3[105-109] | CDR3[114-118] — four distinct blocks, never one continuous strip).

## Context residues are not falsely labeled important

Only true top-5 positions receive the feature-colored highlight band (in both the sequence logo and the quantitative strip below it) and a value printed in the strip. Context positions (the ±2 neighbors) render with plain, unhighlighted logo letters and, where a value happens to exist (i.e. the context position is itself part of the 28-entry empirical interface set), a pale/muted version of the same color scale — never the full-intensity important-position styling.

## VHH/VH logos use identical frequency scales within each feature

Every VHH/VH logo pair within a module is built from the same `pos_label` list, same y-axis limits (0-1), and the same established amino-acid color palette (`figure1_sequence_logo_color_mapping.csv`) — directly comparable by eye, matching the Figure 1/2/9 convention.

## Metric sign conventions (re-verified from source tables, not assumed)

| Feature | Source file | Source column | Verified convention |
|---|---|---|---|
| Sequence divergence | `stage4_position_chemistry_comparison.csv` | `JSD` | unsigned (JSD >= 0 always) |
| Hydropathy shift | `stage4_position_chemistry_comparison.csv` | `delta_mean_KD` | `mean_KD_VHH - mean_KD_VH` (direct recomputation confirms exact match) |
| Charge shift | `Figure9_position_level_bridge.csv` | `delta_charge` | `VHH_mean_sequence_charge - VH_mean_sequence_charge` (confirmed) |
| Native exposure | derived, `Figure9_position_level_bridge.csv` | `VHH_native_exposure - VH_native_exposure` | positive = VHH more exposed |
| CDR3 shielding | `stage6_position_shielding_summary.csv` | `mean_delta_sasa_shield` | VHH-only, no signed comparison attempted |
| Antigen engagement | `stage9_framework_antigen_contact_frequency.csv` | `difference_VHHminusVH` | existing frozen column, confirmed `VHH_antigen_contact_frequency - VH_antigen_contact_frequency` |

## CDR3 shielding is explicitly VHH-only

Panel E's quantitative strip shows only VHH values, uses a distinct sequential violet color scale (never a diverging VHH-vs-VH scale), and the module title is explicitly labeled "CDR3 shielding (VHH only)". No conventional-VH shielding value is shown or implied anywhere in the panel, because no such value exists in the validated source data.

## Unsupported values are NA, not invented

Within the empirical 28-position interface set, a feature value can be genuinely undefined at a given position (e.g. CDR3 shielding is structurally undefined for positions inside CDR3 itself, since CDR3 cannot geometrically shield itself). These cells are rendered as an explicit gray hatched "N/A" cell — confirmed present at position 116 in the CDR3-shielding module (Panel E) and position 116/117 context cells in the antigen-engagement module (Panel F, where the framework-scoped engagement table has no CDR3-range entries). Positions outside the 28-position empirical interface set (most ±2 context positions) simply have no corresponding row in any frozen per-position table and are left blank (thin gray border, no fill, no hatching) — this is a different, weaker claim than "undefined despite being in scope," so it is deliberately rendered differently from the hatched NA cells.

## No heavy scientific calculation was rerun

This rebuild performs only: (1) reading the already-frozen importance ranking from the prior Figure 10 build, (2) reading already-frozen per-position feature values (JSD, delta_hydropathy, delta_charge, delta_native_exposure, CDR3 shielding, antigen-engagement difference) from their original frozen source tables, (3) reading already-frozen cohort amino-acid frequencies from Stage 4, and (4) compositing all of this into six paired sequence-logo + quantitative-strip modules. No interface derivation, sequence curation, FreeSASA, Arpeggio, APBS, PEP-Patch/TPatch, or antigen-contact calculation was performed. No solubility, aggregation, stability, or affinity value was calculated.

## Scope confirmation

Figure 8, Figure 9, and Figure 11 were not modified in this rebuild (verified: file modification timestamps for `Figure8.png`, `Figure9.png`, `Figure11.png` all predate this session's Figure 10 rebuild).
