# Stage 8 — Standard TAP/TNP Author-Code Reproduction Attempt

Per Section 37: attempted to execute the actual TNP author code (`Hydrophobicity_and_Charge_Assigner.CreateAnnotation`) directly on real project structures, to compare against this project's independent reimplementation of the same formula.

## Setup

`external/TNP/lib/python` imports cleanly (package structure intact, no missing Python dependencies beyond numpy/scipy, both present). The real `psa_mac` binary was wired in directly (bypassing the `psa`/`psa_mac` OS-dependent dispatch in `runPSA`, which otherwise expects a `psa` command on `PATH`).

## Result: BLOCKED by a real author-code defect, not a dependency gap

Two independent real VHH structures were tried (`pdb_00001zmy-A`, insertion-coded; `pdb_000028jf-F`, no insertions) through `CreateAnnotation(0, 7.4, <pdb>, 'IG', 'imgt', verbose=False)`. **Both failed identically**:

```
File ".../Hydrophobicity_and_Charge_Assigner.py", line 274, in CreateAnnotation
    res_asa = asa[chain][str(elem[0])+elem[1]]
IndexError: string index out of range
```

This traces to `external/TNP/lib/python/theraprofnano/Hydrophobicity_and_Charge_Profiler/Common/PDBUtils.py`'s residue-key construction in `PDBchain.__init__` (the `sid` splitting logic), which fails on the specific residue-numbering edge cases present in this project's own single-Hchain-only, variable-domain-restricted PDB input — reproducible across both insertion-present and insertion-absent real structures, so the failure is **not** specific to insertion codes (ruling out the most obvious suspect). Root cause not further pursued given it lies inside third-party code not intended to be modified or patched for this project (patching the author's own logic to make it pass would defeat the purpose of an independent validation).

**This is a genuine, documented dependency/compatibility blocker — not a fabricated or approximated agreement.** Per Section 37: "If direct reproduction is impossible due dependency issues: document why. Do not fake agreement." Direct author-code output-level validation is therefore **NOT COMPLETED**.

## What this does and does not affect

- This project's PSH/PPC/PNC implementation (`stage8_patch_scoring.py`) was built by directly reading and transcribing the same source lines of `Hydrophobicity_and_Charge_Assigner.py` (formula, thresholds, normalization, salt-bridge rule) — this is **source-level fidelity**, confirmed and documented in `stage8_tap_source_gate.md`, and is independent of whether the author's own code can be executed end-to-end on this project's specific input format.
- It also independently uses the **same real `psa_mac` binary** as the author code does for the underlying relative-ASA percentage — so the one number that genuinely could not be reconstructed from text (the ASA reference table) is not approximated in this project either; both this project and the (non-executing) author code would, in principle, draw on the same real tool for that value.
- What is **not** obtained is an output-level cross-check (author code's own final PSH/PPC/PNC number vs. this project's number on the same structure) — this specific validation step is honestly reported as blocked, not silently skipped or faked.

**Status: Standard TAP/TNP author-code reproduction — BLOCKED (documented real code defect in third-party PDBUtils.py, not a fabricated pass).**
