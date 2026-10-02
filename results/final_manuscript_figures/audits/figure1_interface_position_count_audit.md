# Figure 1 Revision — 25 vs 28 Interface Position Audit

## A. Authoritative Stage-3 definition

Source: [`src/qc/stage3_position_aggregation.py`](../../src/qc/stage3_position_aggregation.py), thresholds from `configs/stage3.yaml`. A `(imgt_position, insertion_code)` pair is a primary empirical interface member iff:

```
interface_frequency_5  >= 0.25   (ΔSASA >= 5 Å², VH-alone vs VH-with-VL, in >=25% of observations)
observation_fraction   >= 0.25   (observed in >=25% of the 1429 primary-representative structures)
```

Aggregation is grouped by `(imgt_position, insertion_code)` — insertion-code sub-positions are counted as distinct entries, never collapsed into their base integer position.

## B. Exact count in the frozen output

[`data/processed/stage3_empirical_vh_vl_interface_primary.csv`](../../../data/processed/stage3_empirical_vh_vl_interface_primary.csv) (214 total observed `(imgt_position, insertion_code)` rows) contains **exactly 28 rows with `primary_binary_interface == True`**. Verified directly (not from memory) in this session:

```python
sorted 28 (imgt_position, insertion_code) pairs =
[(40,''),(42,''),(44,''),(49,''),(50,''),(52,''),(55,''),
 (66,''),(67,''),(68,''),(69,''),(103,''),
 (107,''),(108,''),(109,''),(110,''),(111,''),(111,'A'),
 (112,''),(112,'A'),(112,'B'),
 (113,''),(114,''),(115,''),(116,''),
 (118,''),(119,''),(120,'')]
```

**25 unique integer IMGT positions**: {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}. This is the exact set literal used throughout the prior "FULL_EMPIRICAL_INTERFACE_28" Python code (`{40,42,44,...,120}`, 25 elements) — that variable name is accurate in describing "28 interface entries" but the set literal itself only ever held the 25 unique integers; the 3 extra entries were always the insertion-coded sub-positions layered on top of positions 111 and 112.

**The 3-entry gap = 111A, 112A, 112B** — three CDR3-loop insertion sub-positions that are independently flagged as interface members (their own ΔSASA/frequency/observation-fraction values cross the same thresholds as their base position), on top of already having 111 and 112 themselves in the set.

## C. Why 3 positions are "missing" from prior sequence-level analyses — root cause

**This is a plotting/filtering omission, not a data-availability or coverage problem.** Checked directly against [`data/processed/stage4_position_chemistry_comparison.csv`](../../../data/processed/stage4_position_chemistry_comparison.csv): **all 28** of the frozen Stage-3 interface entries — including all 3 insertion entries — have real, non-NaN `JSD` and `delta_mean_KD` values (verified: 0 NaN among the 28). Coverage at the 3 insertion entries is lower but well within a usable range, not near-zero:

| position | VHH coverage | VH coverage | JSD | delta_mean_KD |
|---|---|---|---|---|
| 111A | 64.0% (341/533) | 41.4% (591/1429) | 0.0265 | 0.499 |
| 112A | 73.9% (394/533) | 52.0% (743/1429) | 0.0761 | -0.207 |
| 112B | 57.6% (307/533) | 31.8% (455/1429) | 0.0709 | -0.046 |

The earlier Figure 1/2 construction scripts filtered the Stage-4 table with `insertion_code.isna() | (insertion_code=='')` before restricting to the 28-position interface set — a filter meant to collapse deep, mostly-empty CDR3 insertion depths (111C onward, 112C onward, which genuinely are >90% NaN because one or both groups have zero sequences that long) but which, applied indiscriminately, also silently dropped 111A/112A/112B even though those three ARE part of the frozen interface definition and ARE fully sequence-resolved. **This is exactly scenario "accidental plotting omission" from the possible-reasons list — not low coverage, not an upstream definition mismatch, and not a genuine >=50% occupancy filter (Stage 4 applies no such filter at all: `configs/stage4.yaml` sets `coverage_threshold_primary: null` explicitly).**

## D. Definitive statement

**Stage 3 defines 28 empirical interface positions (25 unique IMGT integer positions plus 3 insertion-coded CDR3 sub-positions, 111A/112A/112B, each independently meeting the interface threshold). All 28 are sequence-resolved in the Stage-4 table with real, non-NaN amino-acid-usage and chemistry values — none are excluded for coverage or occupancy reasons. The revised Figure 1 Panel E displays all 28 entries. Figure 2 (main-figure Panel A/B/C, `results/final_manuscript_figures/main/Figure2.png`) currently plots only the 25 base-integer positions due to the same insertion-code filter, and is flagged below as requiring the identical correction — not made in this revision, since the task scope is Figure 1 only.**

## E. Required Figure 2 correction (flagged, not applied)

`results/final_manuscript_figures/main/Figure2.png` (Panels A-C: VHH/VH amino-acid-usage heatmaps and per-position JSD bar chart) was built using the same `insertion_code.isna()` filter and therefore also omits 111A, 112A, and 112B. Per this task's explicit instruction ("Do not automatically modify Figure 2 unless needed for consistency; report the required correction first"), Figure 2 is **not modified in this revision**. The fix, when applied, is mechanical: extend the 28-position axis to include `111A`, `112A`, `112B` as three additional heatmap columns / bar-chart entries between their respective base positions, using the already-existing, already-correct Stage-4 rows (no new computation needed — `stage4_position_aa_frequencies.csv` already has the corresponding `(imgt_position=111, insertion_code='A')` etc. rows populated). This is a one-line filter change in the figure-build script, not a re-analysis.

## F. Full per-position audit table

See [`figure1_empirical_interface_position_audit.csv`](../source_data/figure1_empirical_interface_position_audit.csv) — one row per frozen interface entry (28 rows), with Stage-3 observation/frequency values, Stage-4 sequence-coverage-by-group, prior-Figure-2-inclusion status, and tetrad/CDR3/shieldable-15 annotation flags.
