# Stage 2 Implementation Plan

**Created:** 2026-09-01  
**Scope:** Structure-Level QC and Final Structural Cohort Construction  
**Status:** In Progress

## Architecture Overview

**Modular design:** 1 thin orchestrator + 6 focused helpers + test suite.

```
src/qc/
├── stage2_structure_qc.py       [CREATED] Orchestrator, entry point
├── structure_resolver.py         [CREATED] INSTANCE → mmCIF mapping
├── mmcif_qc.py                  [TODO] mmCIF parsing + entity/chain mapping
├── cdr_qc.py                    [TODO] CDR completeness + IMGT numbering
├── antigen_qc.py                [TODO] Protein-antigen entity identification
├── contact_search.py            [TODO] Distance search + Gordon ≤4.5Å criterion
├── pairing_qc.py                [TODO] VH-VL physical-proximity validation
└── __init__.py

tests/
├── test_structure_resolver.py   [TODO] 6 test cases
├── test_mmcif_qc.py             [TODO] 8 test cases
├── test_contact_search.py       [TODO] 6 test cases (including brute-force validation)
├── test_cdr_qc.py               [TODO] 4 test cases
└── fixtures/
    └── stage2_test_structures/  [TODO] 20 representative mmCIF files for validation
```

## Implementation Priorities (Token-Constrained)

### Phase 1: Contact-Search (CRITICAL)
**Reason:** Computationally most expensive; needs validation before full-cohort run.

**Deliverables:**
- `contact_search.py`: brute-force + optimized (cKDTree) implementations
- `test_contact_search.py`: comparison test on 20-structure set
- Validation report: both algorithms agree within tolerance

**Test set:** 10 VHH–antigen + 10 paired structures with known Gordon contacts

**Expected complexity:** ~500 lines, 2-3 hours runtime for brute-force on 4550 structures

### Phase 2: Helper Modules (In Parallel)
Build remaining helpers in order of dependency:

1. **mmcif_qc.py** — Must precede contact-search (provides structure parsing)
2. **cdr_qc.py** — Depends on mmcif_qc (identifies CDR positions)
3. **antigen_qc.py** — Independent; identifies non-antibody chains
4. **pairing_qc.py** — Depends on mmcif_qc (VH-VL distance)

### Phase 3: Orchestration
- Integrate all helpers into `stage2_structure_qc.py`
- Add checkpointing every 100 structures
- Process full 4550-structure cohort

### Phase 4: Output Generation
- `stage2_structure_qc_master.parquet` (all rows + QC flags)
- Tier tables: `stage2_vhh_qc_{all,strict_primary,secondary_structural}.csv`
- Failure reasons: `stage2_failure_reasons.csv`
- Audit markdown: `stage2_structure_qc.md`

### Phase 5: Validation & Representatives
- Generate 3D interactive HTML (py3Dmol) for 1 VHH–antigen + 1 paired complex
- Create `stage2_representative_structures.csv` with selection justification
- VISUALIZATION DECISION: required (3D geometry interpretation)

## Key Definitions (from configs/stage2.yaml)

| Concept | Definition |
|---------|-----------|
| **CDR completeness (VHH)** | H-CDR1 (27–38), H-CDR2 (56–65), H-CDR3 (105–117) all have ≥1 backbone atom |
| **CDR completeness (paired)** | All 6 CDRs (H1, H2, H3, L1, L2, L3) have ≥1 backbone atom |
| **Backbone atoms** | N, CA, C, O (chemical necessity) |
| **Altloc policy** | Highest occupancy; alphabetical tie-break |
| **Gordon contact** | CDR residue Cα within ≤4.5 Å of protein-antigen Cα |
| **Protein antigen** | Amino-acid polymer (polypeptide(L) or polypeptide(D)); exclude water/ions/heterogens/nucleic acids |
| **VH-VL pairing** | At least one VH atom within 5.0 Å of one VL atom |

## Test-Set Selection

**20 representative structures for brute-force-validation:**

Manually selected from Stage-1 cohort to cover:
- VHH with clear antigen contact (≥1 CDR with contact)
- VHH without clear antigen contact
- Paired VH-VL with both chains present
- Paired VH-VL with only H-chain (misclassified or contamination)
- Structures with missing CDR residues
- Structures with missing CDR1/CDR2/CDR3 (partial VH domain)

**Format:** Extract mmCIF, store in `tests/fixtures/stage2_test_structures/`

## Checkpointing Strategy

Save every 100 structures to `data/interim/stage2/checkpoints/checkpoint_00100.parquet`, etc.

On resume, load latest checkpoint and continue from next structure.

## Output Structure

```
data/interim/stage2/
├── checkpoints/
│   ├── checkpoint_00100.parquet
│   ├── checkpoint_00200.parquet
│   └── ...
└── [final outputs created after full run completes]

results/stage2/
├── stage2_structure_qc_master.parquet    [all 4550 rows + QC flags]
├── stage2_vhh_qc_all.csv                 [all VHH input rows]
├── stage2_vhh_qc_strict_primary.csv      [1328 strict candidates passing QC]
├── stage2_vhh_qc_secondary_structural.csv [mixed/engineered/uncertain passing QC]
├── stage2_paired_qc_all.csv              [all paired input rows]
├── stage2_paired_qc_primary.csv          [3092 conventional pairs passing QC]
├── stage2_failure_reasons.csv            [row, failure_reason_1, reason_2, ...]
├── stage2_failed_structures.csv          [structures that failed QC]
├── stage2_antigen_chain_map.csv          [pdb, chains, antigens, polymers]
├── stage2_missing_residues.csv           [position-level coverage]
├── stage2_structural_coverage_summary.csv [cohort-level stats]
├── stage2_representative_structures.csv   [selection criteria + paths]
├── stage2_structure_qc.md                [audit + methods + findings]
└── structures/
    ├── interactive/
    │   ├── vhh_example_<pdb>_<instance>.html
    │   └── paired_example_<pdb>_<instance>.html
    └── static/
        ├── vhh_example_<pdb>_<instance>.png
        └── paired_example_<pdb>_<instance>.png

results/figures/stage2/
├── stage2_cohort_reduction.png    [IF funnel materially significant]
├── stage2_cohort_reduction.svg
└── figure_manifest.csv
```

## Success Criteria

### Code
- [ ] All 7 modules (orchestrator + 6 helpers) implemented
- [ ] 24 unit tests pass (100%)
- [ ] Full 4550-structure run completes with checkpointing
- [ ] No unhandled exceptions or missing structures

### Data Quality
- [ ] VHH pass rate (% structures passing CDR1/2/3 completeness): ≥85%
- [ ] Paired pass rate (% structures with both H and L chains): ≥90%
- [ ] Gordon contact rate (% structures with ≥1 CDR within 4.5Å of antigen): ≥70%
- [ ] Final VHH structural-QC cohort: 1000–1400 rows
- [ ] Final paired structural-QC cohort: 2500–3000 rows

### Output Validation
- [ ] All tier tables sum to input cohort sizes (no double-counting, no missing rows)
- [ ] Failure reasons justify all exclusions
- [ ] Representative 3D structures render without error
- [ ] Audit markdown is complete and readable

### Documentation
- [ ] Audit markdown explains methods, thresholds, and key findings
- [ ] All figures include sample sizes and semantic colors
- [ ] VISUALIZATION DECISION stated (REQUIRED or NOT REQUIRED with justification)

## Timeline Estimate

- Contact-search + test: **2–3 hours**
- Remaining helpers + test: **4–5 hours**
- Full orchestration + 4550-structure run: **1–2 hours**
- Output generation + audit: **1 hour**
- **Total: 8–11 hours** (wall time; parallel test development)

---

## Next Immediate Steps

1. ✓ Create stage2_structure_qc.py orchestrator
2. ✓ Create structure_resolver.py (INSTANCE → mmCIF)
3. **→ Create contact_search.py (distance algorithm + validation)**
4. Create mmcif_qc.py (mmCIF parsing)
5. Create test suite for all helpers
6. Run contact-search validation on 20-structure test set
7. Integrate all helpers into orchestrator
8. Process full 4550-structure cohort
9. Generate outputs and audit

---

**Status:** Ready to implement Phase 1 (contact-search).

