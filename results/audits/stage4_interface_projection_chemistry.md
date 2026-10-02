# Stage 4 — Empirical Interface Projection and Intrinsic Sequence Chemistry

**Generated:** 2026-09-02
**Scope:** Sequence-level chemistry only (H_intrinsic). No SASA, no CDR3 geometry, no shielding, no mechanistic claims.

---

## A. Strict VHH Structural Input N
1217 (`data/processed/stage2_vhh_strict_primary.parquet`, verified from file, Stage-2 frozen).

## B. Primary Nonredundant VHH Sequence N
**533** redundancy groups from 1217 instances (SABDAB_ID sequence-consistency verified: 300 multi-instance groups checked, 0 inconsistent).

## C. Conventional VH Reference N
**1429** — reused directly from Stage-3's validated paired representative set (`stage3_paired_representatives.parquet`), not re-derived.

## D. VHH Secondary/Sensitivity N
126 (`stage2_vhh_secondary_structural.parquet`), kept fully separate from the primary comparison throughout.

## E. Sequence-Source Validation Result
`_entity_poly_seq` (real-residue truth) cross-referenced with `_atom_site` (resolution status) — the same validated Stage-2/3 mechanism, unmodified. 100% agreement (20/20 spot-checked representatives, all observed residues) between direct atom-level lookup and the `SequenceSource` projection API (`results/audits/stage4_sequence_source_validation.csv`/`.md`).

## F. Projection-Validation Result
**PASS** (`results/audits/stage4_projection_validation.csv`, 10 VHH + 10 VH examples × 214 positions = 2140 checks). All three required states demonstrated (`PRESENT_OBSERVED`, `PRESENT_UNRESOLVED`, `ABSENT_CONFIRMED`), plus the honestly-tracked fourth state `AMBIGUOUS_UNRESOLVED` for cases where IMGT's canonical gap-placement rule would be needed and was not re-derived (that would be equivalent to re-running ANARCI, out of scope). One real bug found and fixed during this validation: N-terminal unresolved stretches (no observed anchor below position 1) were incorrectly falling into `AMBIGUOUS_UNRESOLVED` — fixed by recognizing `_entity_poly_seq`'s `label_seq_id=1` as an always-safe left boundary (nothing can precede the first real residue of a chain).

## G. Primary-Interface Position Coverage Distributions
`results/tables/stage4_interface_projection_coverage.csv`. Zero sequences (VHH or VH) with `primary_interface_position_coverage < 0.5` — the 28-position primary interface is well-covered across the cohort.

## H. Coordinate-Coverage Distributions (Separate)
Tracked in `stage4_sequence_interface_chemistry.{parquet,csv}` as `primary_interface_coordinate_coverage_given_present`, kept strictly distinct from position coverage per spec — a short CDR3 lacking an insertion slot (`ABSENT_CONFIRMED`) never conflated with a present-but-unresolved residue.

## I. H_intrinsic_binary_mean_KD: VHH vs VH
| Group | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| VHH | 533 | **−1.063** | −1.079 | 0.363 | −2.063 | −0.026 |
| VH | 1429 | **−0.696** | −0.696 | 0.371 | −1.629 | 0.861 |

VHH is intrinsically less hydrophobic (lower/more-negative Kyte-Doolittle) than conventional VH at the empirical binary interface, by sequence chemistry alone.

## J. H_intrinsic_weighted_freq5: VHH vs VH
VHH mean **−1.239**, VH mean **−0.703** (full primary sets) — same direction and comparable magnitude to the binary descriptor, confirming the result is not an artifact of the binary threshold choice.

## K. Acidic/Basic/Aromatic/Positive-KD Differences
At the 28-position interface aggregate: `delta_acidic_fraction` mean +0.034, `delta_basic_fraction` mean +0.053 (VHH − VH; both computed per-position then averaged, `stage4_position_chemistry_comparison.csv`). Full per-position breakdown available in that table; no single summary number should be read as more than descriptive.

## L. Nominal Charge Proxy Comparison
`nominal_sidechain_charge = (K+R) − (D+E)` computed per sequence over the binary interface (`stage4_sequence_interface_chemistry.parquet`). Explicitly documented as a formal-charge proxy, not pH-resolved electrostatics (APBS is Stage 8); histidine tracked separately (`histidine_fraction`), never assigned a fixed +1.

## M. Top Interface Positions by JSD
| IMGT pos | Region | JSD | Δmean KD |
|---|---|---|---|
| 50 | FR2 | 0.814 | −7.51 |
| 42 | FR2 | 0.799 | −2.55 |
| 49 | FR2 | 0.721 | −2.35 |
| 52 | FR2 | 0.695 | +2.58 |
| 40 | FR2 | 0.572 | +1.76 |

## N. Top Positions by |Δmean KD|
Position 50 (−7.51), 115/CDR3 (−2.95), 52 (+2.58), 42 (−2.55), 49 (−2.35) — the tetrad plus one CDR3 position (115) and one additional FR2 position (40) dominate both rankings.

## O. Tetrad Residue Distributions and Chemistry Shifts
| Pos | freq≥5 | mean KD VHH | mean KD VH | Δmean KD | JSD | In 28-position interface |
|---|---|---|---|---|---|---|
| 42 | 0.634 | 1.64 | 4.19 | −2.55 | 0.799 | True |
| 49 | 0.965 | −3.06 | −0.72 | −2.35 | 0.721 | True |
| 50 | 1.000 | −3.86 | 3.65 | −7.51 | 0.814 | True |
| 52 | 0.999 | 1.73 | −0.85 | 2.58 | 0.695 | True |

No canonical substitution pattern (e.g. "50=R") was assumed or forced; these are the cohort's actual observed distributions (`stage4_tetrad_comparison.csv`, `stage4_position_aa_frequencies.csv`).

## P. Tetrad-Only vs Full-Interface Result
All 4 tetrad positions fall within the 28-position primary interface, but they carry only **18.5%** of the total `interface_frequency_5` weight and represent 4/28 (14%) of positions — the tetrad is a real subset, not a substitute for the complete empirically-derived interface (`stage4_tetrad_vs_full_interface.csv`).

## Q. Empirical-Interface vs Non-Interface-Framework Descriptive Comparison
Interface (n=28): mean |ΔKD| = 0.958, mean JSD = 0.210. Control framework (n=76, `observation_fraction ≥ 0.25`, CDR excluded by design): mean |ΔKD| = 0.529, mean JSD = 0.160. **This descriptive gap survives the coverage-bias check performed in the mapping audit below (§ mapping bias results)** — reported as descriptive evidence only, no significance test performed (Stage 11).

## R. Primary Representatives vs All-Instance Sensitivity
`results/tables/stage4_primary_vs_all_instance_sensitivity.csv`:

| Group | H_intrinsic_binary_mean_KD (primary) | (all-instance) |
|---|---|---|
| VHH | −1.0625 | −1.0776 |
| VH | −0.6964 | −0.7023 |

Negligible shift (<0.02 in both groups) — repeated crystal structures do not materially alter the primary result.

## S. Secondary VHH Cohort Sensitivity
126 secondary/sensitivity VHH rows (mixed-provenance/engineered/manual-review): mean `H_intrinsic_binary_mean_KD` = **−1.621** (more negative than strict-primary VHH's −1.063). Kept fully separate; no biological claim made about engineered variants shifting toward VH chemistry beyond noting the direction descriptively (`stage4_vhh_secondary_sensitivity.csv`).

## T. pytest Result
**119/119 passed** (full repository suite, including 30 new Stage-4-specific tests).

## U. Generated Figure Paths and Roles
- MAIN_CANDIDATE: `results/figures/stage4/stage4_chemistry_shift.png`/`.svg`
- SUPPLEMENTARY_SEQUENCE_DIVERGENCE: `results/figures/stage4/stage4_positional_divergence.png`/`.svg`
- OPTIONAL_DISTRIBUTION: `results/figures/stage4/stage4_hintrinsic_distribution.png`/`.svg` (generated — distributions cleanly separated with comparable spread, materially aiding interpretation)
- Manifest/README: `results/figures/stage4/figure_manifest.csv`, `results/figures/stage4/README.md`

## V. Limitations
1. H_intrinsic is sequence chemistry only — not actual solvent exposure (H_unshielded/H_native are Stage 5+).
2. A hydrophobic residue by KD scale may be structurally buried; a hydrophilic one may not be surface-exposed.
3. Charge metric is a formal-charge sequence proxy, not pH-resolved electrostatics.
4. No inferential/significance testing performed (Stage 11).
5. The tetrad is one subset of a larger empirically-derived 28-position interface and must not be read as sufficient on its own.
6. ~14 noncanonical (`UNK`) residue observations occur within the primary interface, all from one VH structure (`pdb_00008vxj-A-B`) with a genuinely crystallographer-unassigned CDR3 core — correctly excluded from KD-based statistics, never assigned an invented hydropathy value.

---

## MAPPING-MISSINGNESS BIAS AUDIT

### 1. VHH Ambiguity Rate Within the 28 Primary Interface Positions
**1.74%** (`n_ambiguous_unresolved` / `n_total`, summed across all 28 positions and all 533 VHH representatives).

### 2. VH Ambiguity Rate Within the Same 28 Positions
**2.95%** — modestly higher than VHH, consistent with VH's generally shorter/more-variable CDR3 insertion depth relative to VHH's typically longer CDR3.

### 3. Positions with |Δconfident-mapping-fraction| ≥ 0.05
**5 positions**, all three insertion-coded interface positions plus their immediate insertion-depth neighbors: 111 (Δ=+0.157), 111A (Δ=+0.226), 112 (Δ=+0.104), 112A (Δ=+0.219), 112B (Δ=+0.258) — VHH consistently shows *higher* confident-mapping fraction at these deep CDR3-insertion positions than VH, consistent with VHH's generally longer CDR3 loops making these positions more often genuinely present and resolved.

### 4. Positions with ≥ 0.10 Imbalance
Same 5 positions as above (all 5 exceed 0.10; full detail in `results/audits/stage4_interface_mapping_status_by_position.csv`).

### 5. High-Coverage (≥0.90) Sensitivity Results
| Group | n full | n high-coverage | H_intrinsic_binary_mean_KD (full → high-cov) | H_intrinsic_weighted_freq5 (full → high-cov) |
|---|---|---|---|---|
| VHH | 533 | 394 (74%) | −1.063 → **−1.031** | −1.239 → **−1.223** |
| VH | 1429 | 743 (52%) | −0.696 → **−0.679** | −0.703 → **−0.692** |

Both groups shift only slightly (≤0.03) toward each other under a much stricter, more coverage-matched subset — **the core VHH-vs-VH divergence is not an artifact of unequal mapping confidence.** Primary cohort was not changed based on this threshold (sensitivity only, per instruction).

### 6. Noncanonical Residues Inside the Primary Interface
**14 observations**, all from a single VH structure (`pdb_00008vxj-A-B`), spanning positions 55 and 107–116 (+ insertion variants). Verified directly against the mmCIF: the depositor recorded these residues as `UNK` (unknown/unassignable side-chain density) — a genuine crystallographic annotation, not a parsing artifact. All 14 flagged, none assigned an invented Kyte-Doolittle value; correctly excluded from KD-based statistics and tracked explicitly in `stage4_position_aa_frequencies.csv` under the `X` category.

### 7. Does the Interface-vs-Framework-Control Result Survive Coverage Control?

`results/tables/stage4_interface_vs_control_coverage_check.csv`:

| Category | Group | Mean confident-mapping fraction |
|---|---|---|
| EMPIRICAL_INTERFACE | VHH | 0.937 |
| EMPIRICAL_INTERFACE | VH | 0.900 |
| CONTROL_NONINTERFACE_FRAMEWORK | VHH | 0.994 |
| CONTROL_NONINTERFACE_FRAMEWORK | VH | 0.999 |

Both position sets show VHH ≥ VH coverage (interface gap +0.037; control gap −0.005 — negligible), and both sets are well above 90% confident mapping on average. The interface set is modestly less completely mapped than the framework control in *both* groups (driven by the 3 insertion-coded positions), but the VHH-vs-VH coverage gap itself is small and does not track the much larger chemistry-divergence gap (interface |ΔKD| 0.958 vs control 0.529 — nearly 2×). **Answer: yes, the descriptive result survives the coverage check** — it is not explained by differential VHH/VH mapping confidence between the two position sets.

### 8. Final Primary VHH Representative N
**533**

### 9. Final VH Representative N
**1429**

### 10. Full pytest Result
**119/119 passed**

---

## STAGE 4 STATUS: PASS

**VISUALIZATION DECISION: REQUIRED** — satisfied (1 MAIN_CANDIDATE + 1 SUPPLEMENTARY_SEQUENCE_DIVERGENCE + 1 OPTIONAL_DISTRIBUTION figure, all generated, validated, and manifested).
