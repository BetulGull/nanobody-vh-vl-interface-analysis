# Stage 3 Figures

## FIG_S3_01 — Empirical VH-VL Interface Profile (MAIN_CANDIDATE)

**What it shows:** For every observed VH IMGT position (primary nonredundant
structural cohort), the fraction of structures in which that position had
ΔSASA ≥ 5 Å² due to VL presence (`interface_frequency_5`). Region shading
(FR/CDR), a dashed line at the primary consensus cutoff (0.25), and small
markers for the 4 VHH tetrad positions (42/49/50/52, annotation only — not
treated as predefined interface) are shown. A secondary track below shows
`observation_fraction` per position, since a position observed in very few
structures can otherwise look misleadingly "confident."

**What it does NOT show or prove:** binding energy, thermodynamic stability,
hydrophobicity, or anything about VHH. It is a purely geometric,
population-level empirical summary of 214 VH IMGT positions across
the primary nonredundant paired-structure cohort. 28
positions meet the primary binary consensus definition; this is a starting
definition, not a final biological claim — see the threshold-sensitivity
figure and tables for how it shifts under other reasonable thresholds.

## FIG_S3_02 — Threshold Sensitivity (SUPPLEMENTARY_SENSITIVITY)

**What it shows:** How the number of consensus interface positions changes
across the full grid of per-structure ΔSASA thresholds (1/5/10 Å²),
consensus frequency thresholds (0.10/0.25/0.50), and minimum
observation-fraction thresholds (0.10/0.25/0.50). The primary definition is
one specific cell in this grid, defined *before* looking at these results,
not chosen after inspecting them.

**What it does NOT show or prove:** which threshold is "correct" — no single
threshold combination is presented as uniquely biologically valid. This
figure exists purely to make the sensitivity of the definition transparent.
