# Stage 8 — APBS Full-Cohort Decision

## Benchmark result

`results/audits/stage8_apbs_runtime_benchmark.csv` — 20 real structures (10 VHH including the multi-char-chain/insertion-coded `pdb_00006zbp-FFF`, 10 conventional VH), spanning atom counts 855-1011. All 20 succeeded (`PDB2PQR_SUCCESS_APBS_SUCCESS`). Mean total runtime (PDB2PQR + APBS) = **2.20 s/structure**, tightly clustered (1.95-2.50 s), with no dependence on group (VHH vs VH runtimes overlap fully) and only mild dependence on atom count. No memory issues observed (single-structure sequential runs, ~1GB grid arrays each, transient).

## Projected full-cohort cost

533 VHH + 1429 VH = 1962 structures × ~2.2 s ≈ **~72 minutes** sequential. This is well within the available execution environment (no runtime, memory, or disk blocker identified).

## Decision: **A. FULL COHORT** (preferred, and adopted)

No genuine computational/resource blocker was found. Per Section K of the amendment, the full desired scientific cohort is used:

- **VHH: all 533 primary representatives**
- **VH: all 1429 primary representatives**

No predeclared subsampling is applied. `APBS_EVALUABLE_ANALYSIS_SUBSET` is therefore identical to the full primary cohort (with the small number of individual per-structure failures, if any, tracked explicitly and excluded per-structure, not by design).
