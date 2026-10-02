# Figure 4 — FR2_CDR3_Dis Multimodality Validation Gate

FR2_CDR3_Dis remains the continuous primary variable throughout Stage 6/11 and every downstream figure; this gate only decides whether a bimodal-looking KDE (trough near 9.19 Å, first observed during the earlier exemplar-freeze pass) is robust enough to justify low/high structural exemplars keyed to KDE mode locations, or whether the fallback (deterministic quantile exemplars) should be used instead.

## Data

VHH group only (`FR2_CDR3_Dis` is not computed for VH in the frozen Stage-6 output — see `final_figure_data_lineage.csv`). N=533, range 4.37-14.73 Å, mean=7.95 Å, std=1.99 Å.

## 1. Bandwidth stability scan

Gaussian KDE (`scipy.stats.gaussian_kde`) evaluated at 8 bandwidth-scale factors relative to the default Scott's-rule bandwidth (factor=0.2849):

| bw_scale | bandwidth | n_modes | mode locations (Å) | trough (Å) |
|---|---|---|---|---|
| 0.50 | 0.142 | 5 | 4.46, 6.58, 11.02, 11.78, 14.60 | 4.54 (boundary artifact) |
| 0.70 | 0.199 | 3 | 6.68, 11.03, **14.54** (edge artifact) | 9.09 |
| 0.85 | 0.242 | 3 | 6.73, 11.07, **14.48** (edge artifact) | 9.12 |
| 1.00 (default) | 0.285 | 2 | 6.76, 11.07 | 9.19 |
| 1.15 | 0.328 | 2 | 6.79, 11.05 | 9.27 |
| 1.30 | 0.370 | 2 | 6.80, 11.04 | 9.37 |
| 1.50 | 0.427 | 2 | 6.82, 11.00 | 9.53 |
| 2.00 | 0.570 | 2 | 6.87, 10.63 | 10.29 |

**Not fully stable.** A clean 2-mode structure only appears from bw_scale≥1.0 upward (i.e. at or above the default bandwidth); undersmoothed bandwidths (0.5-0.85) produce 3-5 apparent modes, most of which sit at the extreme tails of the distribution (4.4-4.6 Å, 14.5-14.6 Å) where sample density is thin — these are boundary-kernel artifacts, not evidence of additional real subpopulations. Even within the "stable" 2-mode band (bw_scale 1.0-2.0), the trough location itself drifts from 9.19 to 10.29 Å as bandwidth increases, a >1 Å shift — not a fixed, bandwidth-independent feature.

## 2. Formal diagnostic: Hartigan's dip test

`diptest.diptest(x)` (Hartigan & Hartigan 1985 dip statistic, the standard formal unimodality-vs-multimodality test):

**Dip statistic = 0.02187, p-value = 0.0813**

Not significant at the conventional α=0.05 threshold. The dip test cannot reject the null hypothesis of unimodality for this distribution.

## 3. Group-size stability

Not evaluated further — moot, since the mode/trough structure itself is not bandwidth-robust and the formal test is non-significant; a group split at an unstable threshold would carry unstable group sizes by construction.

## 4. No post-hoc optimization

The trough locations and mode locations above are read directly from unmodified KDE output at predeclared bandwidth-scale factors (a symmetric log-ish sweep 0.5-2.0× the default) — none were selected by checking which threshold best separates shielding or any other downstream outcome. No such downstream check was performed.

## 5. Decision

**Multimodality is NOT robust.** Per the gate's explicit fallback rule, Figure 4's low/high exemplars use **deterministic continuous-quantile selection (Q25/Q75 of FR2_CDR3_Dis)**, not KDE mode locations. This supersedes the earlier exemplar-freeze pass (`final_figure_3d_exemplar_selections.csv`, Figure4 panels D/E), which had used the KDE mode locations (6.76/11.07 Å) before this formal gate was run — those are replaced below with Q25/Q75-based selections for the final-manuscript figure set.

Labels remain **"LOW FR2-CDR3 DISTANCE"** / **"HIGH FR2-CDR3 DISTANCE"** throughout — never "bent"/"extended", consistent with the gate's requirement regardless of which selection rule is used, since neither the earlier KDE-based nor the current quantile-based approach independently confirms two true structural classes in this cohort. FR2_CDR3_Dis itself continues to be treated as a continuous variable in every quantitative panel; only the two illustrative 3D exemplars use a discrete low/high split, and it is presented as an illustrative range comparison, not a claimed subpopulation.

## Updated exemplar selection

VHH group, coverage==1.0, resolution present (N=307 eligible, unchanged eligibility criteria from the original freeze):

| panel | quantile | target (Å) | selected INSTANCE | actual FR2_CDR3_Dis | resolution |
|---|---|---|---|---|---|
| Figure4-D (low) | Q25 | see `figure4_fr2_cdr3_exemplars_quantile.csv` | — | — | — |
| Figure4-E (high) | Q75 | see `figure4_fr2_cdr3_exemplars_quantile.csv` | — | — | — |

Exact frozen values in [`../tables/figure4_fr2_cdr3_exemplars_quantile.csv`](../tables/figure4_fr2_cdr3_exemplars_quantile.csv) (written by `src/qc/final_manuscript_figure4_exemplars.py`).
