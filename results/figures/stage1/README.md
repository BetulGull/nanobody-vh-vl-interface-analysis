# Stage 1 Visualization Summary

This directory contains all primary figures from Stage 1 metadata QC and classification work.

## Figure Organization

Figures are organized by scientific purpose:

- **main/**: Primary scientific figures addressing core questions
- **qc/**: Quality control and diagnostic figures (metadata completeness, redundancy)
- **supplementary/**: (reserved for future supplementary analyses)

Each figure is exported as both PNG (300 dpi) and SVG for maximum flexibility.

## Figure List (Scientific Order)

### FIG_S1_01: Dataset and Curation Flow

**Location:** `main/FIG_S1_01_curation_flow.{png,svg}`

**Question:** How do raw structures partition into primary analysis cohorts?

**Content:**
- Left panel: SD-H classification into six categories (natural camelid, engineered, synthetic, human sdVH, uncertain, other)
- Right panel: Stage-1.5 derived VHH eligibility tiers (strict natural, mixed-provenance, engineered/humanized, manual review)
- Shows the progression from 1597 SD-H metadata rows → 1328 strict natural VHH primary candidates

**Key Finding:** 1419 original metadata natural-camelid candidates reduced to 1328 by downgrading 91 rows from mixed-provenance SABDAB_IDs.

---

### FIG_S1_02: SD-H Classification Distribution

**Location:** `main/FIG_S1_02_sdh_classification.{png,svg}`

**Question:** What are the absolute and relative frequencies of each SD-H metadata category?

**Content:**
- Horizontal bar chart of six SD-H classification categories
- Shows: count, percentage of 1597 total
- Natural camelid highlighted in VHH blue (primary); non-primary in grey

**Key Finding:** Natural camelid candidates dominate (1419 / 1597 = 88.8%) the SD-H raw metadata.

---

### FIG_S1_03: Primary VHH Species Distribution

**Location:** `main/FIG_S1_03_species_distribution.{png,svg}`

**Question:** Which camelid species make up the strict primary VHH cohort?

**Content:**
- Grouped bar chart showing instances vs unique SABDAB_IDs per species
- Species: *Camelus dromedarius*, *C. bactrianus*, *Lama glama*, *Vicugna pacos*, other/unresolved
- Breakdown of 1328 strict primary rows across 1054 unique SABDAB_IDs

**Key Finding:** Llama (*L. glama*, 774 instances, 273 unique IDs) and alpaca (*V. pacos*, 474 instances, 225 IDs) dominate the primary cohort.

---

### FIG_S1_04: FV Construct Architecture

**Location:** `main/FIG_S1_04_fv_construct.{png,svg}`

**Question:** Why are only conventional FV constructs included in the primary paired-antibody cohort?

**Content:**
- Bar chart of FV construct types: conventional (294), SCFV (223), DIABODY (17), CODV (1), DVD (1)
- Conventional highlighted in antibody pink (primary); engineered chain-sharing constructs in grey (excluded)
- Clear annotation: "Only CONVENTIONAL (n=294) is primary paired-antibody cohort"

**Key Finding:** Chain-sharing constructs (SCFV, DIABODY, etc.) are mechanistically distinct from paired-antibody biology; excluded from primary cohort.

---

### FIG_S1_05: Final VHH Provenance Tiers

**Location:** `main/FIG_S1_05_provenance_tiers.{png,svg}`

**Question:** How does the Stage-1.5 mixed-provenance policy stratify the VHH candidates for subsequent analysis?

**Content:**
- Bar chart of four derived eligibility tiers
  - Strict natural VHH (1328): full primary, all analyses
  - Mixed-provenance (108): structural QC + sensitivity only
  - Engineered/humanized (3): secondary cohort
  - Manual review (19): uncertain, requires expert review
- Annotation: "91 of the 108 mixed-provenance rows were previously in the natural-camelid primary cohort; they were downgraded, not excluded"

**Key Finding:** Mixed-provenance policy is conservative: downgrade to lower analysis tier rather than delete; preserve all sequences for comparative/control work.

---

### FIG_S1_06: Redundancy and Pseudoreplication

**Location:** `qc/FIG_S1_06_redundancy.{png,svg}`

**Question:** How many times is the same variable-domain sequence (SABDAB_ID) crystallized across different PDB entries?

**Content:**
- Histogram of instances per SABDAB_ID for the 1328 strict primary VHH rows
- Shows: 1054 unique SABDAB_IDs, 1328 total instances
- Statistics overlay: median=1, mean≈1.26, max shown
- Most IDs appear once; a few heavily repeated

**QC Purpose:** Quantify pseudoreplication before Stage 2 structure analysis. For downstream work, consider one representative per SABDAB_ID to avoid non-independence.

---

### FIG_S1_07: Metadata Quality Assessment

**Location:** `qc/FIG_S1_07_metadata_quality.{png,svg}`

**Question:** What is the quality and completeness of crystallographic metadata across cohorts?

**Content:**
- Four-panel figure
  - Panel 1: resolution distribution by cohort (SD-H blue, FAB pink, FV purple)
  - Panel 2: r_factor data completeness (% with values)
  - Panel 3: r_free data completeness (% with values)
  - Panel 4: sample size by cohort (n=1597 SD-H, 2798 FAB, 536 FV)

**QC Purpose:** Document metadata readiness before coordinate-level Stage 2 work. All three cohorts pass basic structural criteria (X-ray, resolution ≤3.0 Å); additional missingness in r_factor/r_free noted.

---

### FIG_S1_08: Mixed-Provenance Audit

**Location:** `main/FIG_S1_08_mixed_provenance_audit.{png,svg}`

**Question:** What is the composition of each mixed-provenance SABDAB_ID?

**Content:**
- Five-panel figure (one per SABDAB_ID)
- Each panel: horizontal bar breakdown showing which Stage-1 classification labels apply to instances within that ID
  - H00QEL0000: 13 instances → natural camelid vs human sdVH
  - H021EL0000: 23 instances → natural camelid vs human sdVH
  - H02DAL0000: 2 instances → natural camelid (VHH81) vs human sdVH (caplacizumab, humanized derivative)
  - H02Y5L0000: 11 instances → natural camelid vs human sdVH (VHH + MHC antigen complexes; Stage-1.5 confirmed no contamination)
  - H06BKL0000: 2 instances → Vicugna pacos both, but heavy_taxid inconsistency (one row self-contradictory text/taxid)

**Key Annotation:** "Actual variable-domain sequences were identical at shared IMGT positions; disagreement is provenance metadata, not sequence difference."

**Audit Conclusion:** Provenance inconsistency is real; affects only these 5 IDs (108 rows) — conservative approach is to downgrade from primary until resolved, but never to exclude or distrust sequence identity.

---

## Summary Statistics

| Metric | Count |
|--------|-------|
| SD-H rows in raw metadata | 1597 |
| SD-H natural camelid candidates | 1419 |
| Strict primary natural VHH (Stage-1.5) | 1328 |
| Mixed-provenance VHH (downgraded) | 108 |
| Unique SABDAB_IDs in strict primary | 1054 |
| Engineered/humanized secondary | 3 |
| Manual review (uncertain) | 19 |
| FAB rows (all conventional paired) | 2798 |
| FV conventional paired | 294 |
| FV total | 536 |
| **Total dataset size** | **4931** |

## Reproducibility

All figures are generated programmatically from:

- **Source:** `src/visualization/stage1_figures.py`
- **Config:** `configs/visualization.yaml` (semantic color palette)
- **Data:** `data/metadata/master_stage1_resolved.parquet`

To regenerate:
```bash
micromamba run -n vhh3d python3 src/visualization/stage1_figures.py
```

## Color Palette

All figures use the project-standard **Gordon-inspired** semantic color scheme:
- VHH / sdAb: **#4D69B2** (blue)
- Conventional antibody: **#EF7F7E** (pink)
- Conventional VH: **#8670B1** (purple)
- Grey: **#CCCCCC** (missing/excluded)

See `configs/visualization.yaml` and `docs/visualization_standard.md` for full specification.

## File Formats

- **PNG:** 300 dpi, white background, suitable for publication
- **SVG:** vector format, editable with Inkscape/Adobe/etc.

## Next Steps

Stage 1.5 concludes with this visualization suite. All figures undergo programmatic validation (dimensions, format, consistency) before `STAGE 1 VISUALIZATION STATUS: PASS`.

Stage 2 (structure-level QC and mechanistic analysis) will require additional figures addressing:
- Coordinate completeness and QC by chain/domain
- Antigen contact coverage
- CDR3 geometry and shielding
- Former VH–VL interface empirical definition
- 3D molecular visualizations (static + interactive)

---

**Generated:** 2026-09-01  
**Status:** Complete
