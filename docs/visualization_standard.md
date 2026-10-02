# Visualization Standard — VHH Structural Mechanistics Project

**Version 1.0.0 | Updated 2026-09-01**

This document specifies the visual standards for all figures generated in this project, from Stage 1 through Stage 2+ mechanistic analysis. The goal is consistency, reproducibility, and clarity in communicating structural and sequence data.

**Configuration source:** `configs/visualization.yaml` (machine-readable; this document is the human narrative.)

## 1. Semantic Color Palette

All figures use a **project-standard semantic color scheme** inspired by Gordon et al. (2023) representations of immunoglobulin superfamily domains. This ensures that the same biological entity is visually recognized across all figures and publications.

### Core Semantic Colors

| Concept | Color | Hex | Usage |
|---------|-------|-----|-------|
| VHH / sdAb | Gordon blue | #4D69B2 | Single-domain antibody, primary representation |
| VHH (secondary) | Light blue | #7D92C9 | 2nd series, background, lower priority |
| Conventional antibody | Gordon pink | #EF7F7E | Fab / full-length antibody, VL |
| Conv. antibody (secondary) | Light pink | #F4A5A8 | 2nd series, lower priority |
| Conventional VH | Gordon purple | #8670B1 | VH from paired VH–VL structures |
| Protein antigen | Restrained green | #6BAA63 | Target / bound antigen |
| CDR3 focus | Cyan / teal | #00BCD4 | Mechanistic highlight, shielding |
| Empirical interface | Yellow-orange | #FFB300 | Former VH–VL contact surface |
| Missing / excluded | Neutral grey | #CCCCCC | Reference, excluded, missing data |

### Quantitative Scales

When representing continuous data, use:

- **Hydrophobic exposure:** yellow → orange → red (sequential scale)
- **Electrostatic potential (APBS):** blue (positive) and red (negative), always with labeled diverging scale
  - **CRITICAL:** electrostatic blue ≠ VHH semantic blue. Always include a color bar or legend to distinguish.

## 2. Figure Standards

### Resolution and Format

Every scientific figure must be exported as:

- **PNG: 300 dpi minimum** (for print/publication)
- **SVG: vector format** (for editing, re-scaling)
- **PDF: vector** (when straightforward; optional for supplementary)

### Layout and Design

- **Background:** white
- **Typography:** sans-serif, readable at 6 points when printed
- **Titles:** bold, descriptive, state the actual biological question
- **Legends:** consistent formatting across all figures; do not cover data
- **Axes:** clear labels; do not truncate without justification and explicit notation

### Data Representation

1. **Always show sample size (N)** in axis labels, legend, or caption
2. **Show distributions, not just means** where relevant (use histograms, violin plots, ECDF, density, hexbin, or scatter with transparency)
3. **Handle category imbalance:**
   - Use log scale if needed
   - Use inset plots or split panels for rare categories (e.g., CODV, DVD)
   - Annotate tiny categories explicitly so they don't disappear
   - Never allow large categories to visually erase small ones
4. **Remove chart junk:** no decorative gradients, unnecessary 3D, or misleading visual tricks

## 3. Figure Hierarchy

Figures are organized by scientific purpose:

```
results/figures/stage1/
├── qc/                # data quality, missingness, coverage
├── main/              # primary scientific figures
└── supplementary/     # supporting analyses

results/figures/stage2/  # future structure-level figures
├── qc/
├── main/
└── supplementary/

results/structures/     # future 3D molecular visualizations
├── static/            # PNG/PDF static molecular images
├── interactive/       # standalone HTML 3D viewers
└── scripts/           # PyMOL/ChimeraX sessions, regeneration scripts
```

## 4. Reproducibility and Provenance

Every figure must be reproducible from code. For each figure, retain:

- **source_script:** Python/R/other that generates the plot
- **input_table:** CSV/Parquet path
- **filters_applied:** describe any row filtering
- **variables_shown:** which columns are plotted
- **sample_size:** N used (must match the audited table)
- **random_seed:** if applicable (e.g., for jitter, subsampling)

**Never:**
- Manually edit scientific values into images
- Hard-code values from analysis into figure generation
- Use different colors for the same concept across different figures
- Cherry-pick "impressive" structures without documented reproducible selection criteria

**Always:**
- Load colors from `configs/visualization.yaml` (never hard-code)
- Programmatically generate annotations from audited metadata tables
- Validate that generated figures actually exist and are readable before declaring pass

## 5. Stage 1: Metadata and QC Figures

Stage 1 produces 8+ required figures documenting data provenance, classification, and quality.

### FIG_S1_01: Dataset and Curation Flow

**Question:** How many structures entered at each stage, and how do they partition into primary/secondary cohorts?

**Content:**
- Downloaded source cohorts: SD-H (1597), FAB (2798), FV (536)
- SD-H classification: natural camelid, engineered/humanized, synthetic, human sdVH, other, uncertain
- Derived Stage-1.5 VHH tiers: strict natural, mixed-provenance, engineered secondary, manual-review, excluded
- Conventional paired antibodies: FAB (all conventional, 2798), FV breakdown (conventional 294, SCFV 223, diabody 17, CODV 1, DVD 1, etc.)
- Visual: Sankey/alluvial (if legible) or multi-panel funnel/flow showing row counts at each step
- Must show: why only conventional FAB and conventional FV enter the primary paired-antibody cohort

### FIG_S1_02: SD-H Classification Distribution

**Question:** What is the breakdown of Stage-1 SD-H classification categories?

**Content:**
- Six original SD-H categories (natural camelid, engineered/humanized, synthetic, human sdVH, other, uncertain)
- Show: count, percentage, sample size
- Use VHH blue prominently for natural camelids
- Use neutral grey for non-primary categories
- Visual: horizontal bar chart or side-by-side bars with N labeled

### FIG_S1_03: Primary VHH Species Distribution

**Question:** Which camelid species dominate the natural-VHH primary cohort?

**Content:**
- Species groups: Camelus dromedarius, C. bactrianus, C. ferus, other Camelus, Lama glama, other Lama, Vicugna pacos, other/unresolved Camelidae
- Show both: instance count AND unique SABDAB_ID count per species
- Visual: paired horizontal bars, dot plot, or back-to-back bar chart
- Make clear which counts are instances vs. unique sequences

### FIG_S1_04: FV Construct Architecture

**Question:** Why does construct type matter for paired-antibody primary cohort?

**Content:**
- FV construct categories: conventional (294), SCFV (223), DIABODY (17), CODV (1), DVD (1), other
- Make visually obvious why only conventional FV is primary (chain-sharing constructs exclude)
- Visual: stacked bar, pie, or grouped bar highlighting the dominating conventional category
- Annotate small categories (CODV, DVD) explicitly

### FIG_S1_05: Final VHH Provenance Tiers

**Question:** How does the Stage-1.5 mixed-provenance policy partition the natural-VHH candidates?

**Content:**
- Derived VHH structural-QC pool decomposition:
  - 1328 strict natural camelid (primary)
  - 108 mixed-provenance same-sequence (downgraded)
  - 3 engineered/humanized secondary
  - 19 manual-review/uncertain
- Critically: show that 91 of the 108 mixed-provenance rows were *previously* inside the natural-camelid primary set
- Visual: stacked bar or grouped bar showing partition, with callout or annotation for the 91 "removed" rows
- Purpose: make clear that mixed-provenance is a downgrade, not an exclusion

### FIG_S1_06: Redundancy and Pseudoreplication

**Question:** Are certain variable-domain sequences (SABDAB_IDs) heavily duplicated?

**Content:**
- Distribution of instances per SABDAB_ID for:
  - VHH-side candidates
  - Conventional paired-antibody candidates
- Visual: ECDF, histogram, or log-scaled distribution
- Label or annotate highest-repetition IDs
- Do not allow a few highly repeated IDs to dominate the visualization (use log scale if needed)

### FIG_S1_07: Metadata Quality

**Question:** What is the structural-metadata quality distribution?

**Content:**
- Resolution distribution by source cohort (SD-H, FAB, FV)
- R-factor missingness / presence
- R-free missingness / presence
- Visual: overlaid histograms or violin plots by cohort; separate panel for R-factor/R-free if present
- Purpose: document metadata completeness before Stage 2 coordinate-level work

### FIG_S1_08: Mixed-Provenance Audit

**Question:** What is the nature of the 5 mixed-provenance SABDAB_IDs?

**Content:**
- For each ID: number of rows, species/taxid labels, prior classification, final derived status
- Clear annotation: "Actual variable-domain sequences were identical at shared IMGT positions; disagreement is provenance metadata, not sequence"
- Do NOT imply sequence discrepancy
- List the 5 IDs and their core finding (e.g., H02Y5L0000 has Lama glama instances + human instances, same sequence)
- Visual: table or annotated multi-panel figure showing each ID's composition

## 6. Future Structural Visualization (Stage 2+)

### Static Molecular Representations

All 3D structures must be visualized as:

- **Static PNG/PDF:** 300 dpi, labeled with PDB, INSTANCE, SABDAB_ID, chains
- **Components clearly shown:** Hchain, Lchain (if present), antigen
- **Stereochemistry:** correct and inspectable

### Interactive 3D Molecular Viewers

For structural stages, generate interactive views:

- **Technology:** py3Dmol / 3Dmol.js standalone HTML
- **Capabilities:** rotation, zoom, residue inspection, chain selection
- **Labeling:** PDB, INSTANCE, SABDAB_ID in title
- **Reproducibility:** preserve session file or regeneration script

### Representative Structure Selection

Structures shown as examples must be chosen reproducibly:

- **Valid selection criteria:** median metric, 25th/75th percentile, high/low specific measure, matched control
- **Invalid:** visual cherry-picking, aesthetic impressiveness
- **Documentation:** results/structures/representative_selection.csv with justification

### Stage 2 3D Examples (when reached)

Require at least:

1. Valid VHH + protein antigen complex
2. Valid paired VH–VL + protein antigen complex
3. Major QC failure classes if present (e.g., missing backbone, ambiguous antigen)

Each labeled clearly: PDB, INSTANCE, SABDAB_ID, chains, cohort.

## 7. Special Visualization Requirements by Analysis

### Former VH–VL Interface (Stage 2+)

When the empirical former-VL-interacting surface is derived:

- **Generate paired interactive views:**
  - Conventional VH + VL (full complex)
  - Same VH with VL hidden
  - Same VH with interface positions highlighted
- **Color interface by:**
  - Interface frequency (how many VH–VL structures contact this position)
  - Mean Δ SASA per position
  - Median Δ SASA per position
- **Later project positions onto VHH** to compare
- **Critical:** interface must come from empirical Δ SASA, NOT assumed to be FR2

### CDR3 Shielding (Stage 2+)

When CDR3 shielding metrics are calculated:

- **Generate paired views:**
  - Native VHH (actual structure)
  - Same coordinates with CDR3 omitted (geometric counterfactual)
- **Label the second view explicitly:** "geometric counterfactual — not a physical folded structure"
- **Highlight:** CDR3, former VH–VL interface projection, newly exposed surface
- **Accompany with:** Δ SASA shielding numbers, shielding fraction, Δ H_CDR3

### CDR3–Interface Packing (Arpeggio, Stage 2+)

- **Distinguish interaction types:** H-bonds, hydrophobic, aromatic, ionic, vdW
- **Do NOT present generic proximal contacts as chemically specific**
- **Use clear glyphs or labels** for each interaction type

### Patch Topology and Electrostatics (Stage 2+)

Maintain **visual separation:**

- **PPC/PNC discrete charge-patch descriptors:** discrete patch locations (color by patch type or magnitude)
- **APBS continuum electrostatic potential:** continuous potential surface (use diverging blue–red scale with labeled scale bar)
- **Never merge** patch and electrostatic into generic "charge" visualization

## 8. Validation Checklist

Before a figure is considered complete:

- [ ] PNG file exists and dimensions are nonzero
- [ ] SVG file exists and is valid XML
- [ ] Source script is documented and reproducible
- [ ] Input table path is correct and file exists
- [ ] Sample size (N) matches the audited source table
- [ ] All labels are readable (not truncated)
- [ ] Title states the biological question
- [ ] Legend does not cover data
- [ ] Categories with small N are not invisible
- [ ] Semantic colors loaded from configs/visualization.yaml
- [ ] No manually edited scientific values
- [ ] Random seed documented (if applicable)

## 9. Figure Manifest and README

Every stage produces:

- **figure_manifest.csv:** machine-readable index with figure_id, filename, title, source_table, n_used, etc.
- **README.md:** human-readable summary showing figures in order with explanatory captions

## 10. Scientific Visualization Gate (Stage 2+)

**Key principle:** Generate a figure only when it serves the science, not because data exists.

From Stage 2 onward, create a figure ONLY when at least ONE of these is true:

1. **Direct mechanistic contribution:** The figure directly addresses the central question (e.g., VHH CDR3 shielding vs. conventional VH burial)
2. **Major QC decision:** A dataset/cohort change is more clearly understood visually than in tabular form
3. **Conventional interpretation:** The quantity is spatially, distributionally, or structurally interpreted in the literature (e.g., molecular surfaces, geometry distributions)
4. **Structural necessity:** A 3D view materially improves interpretation of the result

**Do NOT automatically visualize:**
- Missingness or completeness percentages
- Species counts or demographic breakdowns (unless addressing the central question)
- Software checks or method validation tables
- Every intermediate metric
- Every exclusion reason
- Metadata fields simply because a CSV exists

**No figure quota:** A stage may have 0, 1, or many figures. Fewer is better if tables suffice.

**VISUALIZATION DECISION per stage:** REQUIRED (justified by scientific question) or NOT REQUIRED (table is sufficient)

Fail visualization ONLY when a scientifically necessary figure was identified but not successfully produced.

## 11. Overall Project Visualization Status

Stage 1 receives `VISUALIZATION STATUS: PASS` if:
- All required figures exist and are valid
- Sample sizes are shown
- Distributions are shown where relevant
- Colors are used consistently per config
- Provenance is documented
- Validation checklist is complete

**Stage 2+ receives PASS if:**
- VISUALIZATION DECISION is explicitly stated (REQUIRED or NOT REQUIRED)
- If REQUIRED: the scientific question is stated, and figures exist with that scientific role
- If NOT REQUIRED: tabular/audit output is sufficient and documented
- All generated figures satisfy the scientific visualization gate
- Figure roles (MAIN/SUPPLEMENTARY/QC_ONLY) are assigned

---

## 11. Stage-Wise Visualization Policy (CORRECTED SCOPE)

### Stage 1: Metadata Audit
**VISUALIZATION DECISION:** REQUIRED  
**Figures generated:** 8 (curation, classification, species, constructs, tiers, redundancy, metadata quality, provenance)  
**Role of figures:** Communicate dataset assembly and cohort composition to readers.

| Figure | Role | Rationale |
|--------|------|-----------|
| FIG_S1_01 (Curation Flow) | MAIN_CANDIDATE | Core data-flow visualization; communicates cohort assembly and classification decisions to readers |
| FIG_S1_02 (SD-H Classification) | SUPPLEMENTARY_CURATION | Curation context; detail redundant with FIG_S1_01 |
| FIG_S1_03 (Species Distribution) | SUPPLEMENTARY | Diversity context; not mechanistically central |
| FIG_S1_04 (FV Construct) | SUPPLEMENTARY_QC | QC rationale; explains paired-antibody filtering |
| FIG_S1_05 (Provenance Tiers) | SUPPLEMENTARY | Likely redundant with FIG_S1_01; candidate for consolidation |
| FIG_S1_06 (Redundancy) | QC_ONLY | Diagnostic for deduplication; informs Methods section |
| FIG_S1_07 (Metadata Quality) | QC_ONLY | Structural metadata readiness; likely in Methods or supplement |
| FIG_S1_08 (Mixed-Provenance Audit) | QC_ONLY | Provenance conflict audit; detailed standalone audit or supplement |

### Stage 2: Structure-Level QC and Cohort Validation
**VISUALIZATION DECISION:** CASE-BY-CASE (per scientific gate)  
**Scope:** Hchain/Lchain mapping, variable-domain completeness, IMGT coverage, CDR1/2/3 presence, paired validation, antigen identification, Gordon-compatible antigen-contact filtering, cohort reduction accounting.  
**Figures justified by scientific gate ONLY IF:**
1. Structure-level QC materially reduces the dataset → cohort-reduction funnel (one figure, table-alternative)
2. Spatial interpretation of antigen-binding geometry required → representative 3D complex examples (validated VHH–antigen + validated paired VH–VL–antigen)
3. Recurrent QC-failure structure discovered with mechanistic significance → document failure example(s)

**Figures NOT justified in Stage 2:**
- CDR3 geometry distributions (that is Stage 3+: mechanistic geometry analysis)
- CDR3 shielding plots (that is Stage 3+: mechanistic shielding analysis)
- Empirical VH–VL interface derivation/projection (that is Stage 3+: interface mechanistics)
- Hydrophobic-surface analysis (that is Stage 3+: surface chemistry)
- Arpeggio interaction networks (that is Stage 3+: packing mechanistics)
- APBS electrostatics (that is Stage 4+: electrostatic mechanistics)

### Stage 3+: Mechanistic Analysis
**VISUALIZATION DECISION:** CASE-BY-CASE (per scientific gate)  
**Scope includes:** CDR3 geometry, shielding, interface topology, hydrophobic patches, Arpeggio packing, electrostatics.  
**Visualization justified:** Spatial/distributional/mechanistic interpretation is central to all of these; 3D views likely needed.

---

**For questions or revisions to this standard, refer to:**
- `configs/visualization.yaml` — authoritative settings (includes scientific visualization gate)
- `results/figures/stage1/figure_manifest.csv` — figure metadata and roles
- `results/figures/` — all stage figures and audit outputs
