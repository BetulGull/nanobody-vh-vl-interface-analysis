# Stage 4 Figures

## FIG_S4_01 — Intrinsic Sequence Chemistry Shift (MAIN_CANDIDATE)

**What it shows:** For each of the 28 Stage-3 primary empirical
VH-VL interface positions, the difference in mean Kyte-Doolittle intrinsic
hydropathy between VHH and conventional VH sequences at that exact IMGT
position. Point size encodes the Stage-3 `interface_frequency_5` (how
consistently that position is buried by VL across paired structures); the
four VHH tetrad positions (42/49/50/52) are annotated in teal — as
annotation only, not as a predefined or privileged subset.

**What it does NOT show or prove:** solvent-exposed hydrophobicity (this is
sequence-level Kyte-Doolittle chemistry only — no SASA calculated in this
stage), binding energy, aggregation propensity, or anything about VHH
stability/solubility directly. A negative Δmean KD means VHH is
intrinsically less hydrophobic (by sequence alone) than conventional VH at
that position — nothing more.

## FIG_S4_02 — Positional Divergence (SUPPLEMENTARY_SEQUENCE_DIVERGENCE)

**What it shows:** Jensen-Shannon divergence between VHH and VH amino-acid
distributions at each of the 28 primary interface positions, distinguishing
framework from CDR positions, with tetrad positions marked.

**What it does NOT show or prove:** statistical significance — no
multiple-testing-corrected inference is performed in Stage 4 (deferred to
Stage 11). This is a descriptive ranking only.

## FIG_S4_03 — H_intrinsic Distribution (OPTIONAL_DISTRIBUTION)

**What it shows:** Full distributions (not just means) of
`H_intrinsic_binary_mean_KD` and `H_intrinsic_weighted_freq5` for the
primary nonredundant VHH (n=533) and VH (n=1429) representative sets,
confirming the shift seen in FIG_S4_01 is a population-wide pattern.

**What it does NOT show or prove:** causality, mechanism, or that this
reflects actual solvent exposure in either group.
