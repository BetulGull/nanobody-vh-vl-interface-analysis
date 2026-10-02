"""
Stage 3 — FreeSASA Validation Gate (spec §23-25)

Validates, on a deterministic real-structure set, before any full-cohort run:
  D. per-residue VH SASA sum ≈ total VH SASA from the same calculation
  E. repeated run gives identical values (reproducibility)
  F. CLI vs Python-API comparison: NOT attempted (documented reason below) —
     using internal reproducibility + the geometric translated-VL control
     instead, per the spec's explicit allowance.
  §24: translated-VL geometric control — VL translated far away must leave
     VH SASA ≈ VH-alone SASA (ΔSASA ≈ 0)
  §25: total ΔSASA consistency — sum_i(delta_i) ≈ total_VH_alone - total_VH_with_VL

CLI-vs-API rationale for NOT forcing comparison: the CLI operates on a
fixed-column PDB text format (atom name justification, chain ID single
character, insertion code column) reconstructed from our in-memory atom
list; producing a byte-correct PDB subset (VH-variable-domain-only, with the
exact same altloc/insertion-code selections already applied) introduces its
own formatting-fidelity risk that would need separate validation to trust —
i.e. comparing against a CLI run would validate the PDB writer, not the SASA
calculation itself. The two methods available here (Python API `addAtom`
direct coordinate injection, used identically twice for the reproducibility
check) and the geometric translated-VL control are cleaner, more direct
validations of the actual computation than an indirect CLI round-trip.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
import copy

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa, get_freesasa_config_record

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

RESIDUE_SUM_TOLERANCE = 1e-6  # Å² -- exact within float precision, same underlying calc
REPRODUCIBILITY_TOLERANCE = 1e-9  # Å² -- identical deterministic algorithm, no RNG
TRANSLATED_VL_TOLERANCE = 0.5  # Å² -- allow tiny numerical noise, no real contact expected


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def validate_structure(instance: str, master: pd.DataFrame) -> dict:
    row = master[master['INSTANCE'] == instance].iloc[0]
    hchain = row['Hchain']
    lchain = normalize_lchain(row['Lchain'])
    path = PROJECT_ROOT / row['extracted_structure_path']

    s = MMCIFStructure(str(path))
    model_id = s.get_model()

    vh_atoms = extract_variable_domain_atoms(s, hchain, model_id)
    vl_atoms = extract_variable_domain_atoms(s, lchain, model_id) if lchain else []

    result = {'INSTANCE': instance}

    # --- D: per-residue sum vs total, VH-alone ---
    per_res_alone, total_alone = calc_sasa(vh_atoms)
    sum_per_res_alone = sum(per_res_alone.values())
    result['vh_alone_total_sasa'] = total_alone
    result['vh_alone_sum_per_residue'] = sum_per_res_alone
    result['vh_alone_residue_sum_consistent'] = abs(total_alone - sum_per_res_alone) < RESIDUE_SUM_TOLERANCE

    # --- D: per-residue sum vs total, VH+VL ---
    per_res_with_vl, total_with_vl = calc_sasa(vh_atoms + vl_atoms)
    # total_with_vl includes VL atoms' own SASA too; for VH-only sum we need VH residues only
    vh_chain_resnums = set((a['auth_asym_id'], f"{a['auth_seq_id']}{a['ins_code']}") for a in vh_atoms)
    vh_total_within_complex = sum(v for k, v in per_res_with_vl.items() if k in vh_chain_resnums)
    result['vh_with_vl_sum_per_residue_vh_only'] = vh_total_within_complex

    # --- E: reproducibility (repeat calc) ---
    per_res_alone_2, total_alone_2 = calc_sasa(vh_atoms)
    result['reproducible_total'] = abs(total_alone - total_alone_2) < REPRODUCIBILITY_TOLERANCE
    result['reproducible_per_residue'] = all(
        abs(per_res_alone[k] - per_res_alone_2[k]) < REPRODUCIBILITY_TOLERANCE for k in per_res_alone
    )

    # --- §24: translated-VL geometric control ---
    if lchain and vl_atoms:
        vl_far = copy.deepcopy(vl_atoms)
        for a in vl_far:
            a['x'] += 100000.0  # far enough that no contact is possible; in-memory only, never saved
        per_res_far, total_with_far_vl = calc_sasa(vh_atoms + vl_far)
        vh_total_with_far_vl = sum(v for k, v in per_res_far.items() if k in vh_chain_resnums)
        delta_far = total_alone - vh_total_with_far_vl
        result['translated_vl_vh_sasa'] = vh_total_with_far_vl
        result['translated_vl_delta_from_alone'] = delta_far
        result['translated_vl_control_pass'] = abs(delta_far) < TRANSLATED_VL_TOLERANCE
    else:
        result['translated_vl_vh_sasa'] = None
        result['translated_vl_delta_from_alone'] = None
        result['translated_vl_control_pass'] = None

    # --- §25: total ΔSASA consistency ---
    if lchain:
        sum_delta_i = sum(
            per_res_alone[k] - per_res_with_vl.get(k, 0.0)
            for k in vh_chain_resnums
        )
        total_delta = total_alone - vh_total_within_complex
        result['sum_delta_i'] = sum_delta_i
        result['total_delta_from_totals'] = total_delta
        result['total_delta_consistent'] = abs(sum_delta_i - total_delta) < RESIDUE_SUM_TOLERANCE
    else:
        result['sum_delta_i'] = None
        result['total_delta_from_totals'] = None
        result['total_delta_consistent'] = None

    return result


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")

    validation_instances = [
        'pdb_000011ha-G-K',   # FAB, VH+CH1 same chain
        'pdb_00001a14-H-L',   # paired FV
        'pdb_000013fj-H-L',   # insertion-coded
        'pdb_00002xqb-H-L',   # missing coordinates
        'pdb_00001lk3-H-L',   # different resolution / has MSE elsewhere in file (antigen)
    ]

    results = [validate_structure(inst, master) for inst in validation_instances]
    df = pd.DataFrame(results)

    out_path = AUDITS_DIR / "stage3_sasa_validation.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}")
    logger.info(df.to_string(index=False))

    logger.info(f"\nFreeSASA config: {get_freesasa_config_record()}")

    all_pass = (
        df['vh_alone_residue_sum_consistent'].all() and
        df['reproducible_total'].all() and
        df['reproducible_per_residue'].all() and
        df['translated_vl_control_pass'].dropna().all() and
        df['total_delta_consistent'].dropna().all()
    )
    logger.info(f"\nSASA VALIDATION GATE: {'PASS' if all_pass else 'FAIL'}")
    return all_pass, df


if __name__ == "__main__":
    ok, _ = main()
    sys.exit(0 if ok else 1)
