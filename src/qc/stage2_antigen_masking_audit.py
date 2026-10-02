"""
Stage-2 Closure Audit — Multi-Instance PDB Antigen Masking

Verifies that, after the fix, no antibody chain belonging to ANY SAbDab2
INSTANCE sharing a PDB is ever assigned as a candidate protein antigen for
another INSTANCE in that same PDB.

Output: results/audits/stage2_multiinstance_antigen_masking.csv
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")

    pdb_to_ab_chains = {}
    for _, row in master.iterrows():
        pdb = row['PDB']
        h = row['Hchain']
        l = normalize_lchain(row['Lchain'])
        pdb_to_ab_chains.setdefault(pdb, set()).add(h)
        if l:
            pdb_to_ab_chains[pdb].add(l)

    pdb_instance_counts = master.groupby('PDB')['INSTANCE'].nunique()
    multi_instance_pdbs = pdb_instance_counts[pdb_instance_counts > 1]

    rows = []
    n_false_assignments = 0
    n_instances_affected = 0

    for _, row in master.iterrows():
        pdb = row['PDB']
        if pdb not in multi_instance_pdbs.index:
            continue

        own_h = row['Hchain']
        own_l = normalize_lchain(row['Lchain'])
        own_chains = {own_h}
        if own_l:
            own_chains.add(own_l)

        other_ab_chains = pdb_to_ab_chains.get(pdb, set()) - own_chains
        antigen_chains = [c for c in str(row.get('antigen_chains', '') or '').split(';') if c]

        false_assignment_chains = set(antigen_chains) & other_ab_chains

        if false_assignment_chains:
            n_false_assignments += 1
            n_instances_affected += 1

        rows.append({
            'PDB': pdb,
            'INSTANCE': row['INSTANCE'],
            'own_antibody_chains': ';'.join(sorted(own_chains)),
            'other_instances_antibody_chains_in_same_pdb': ';'.join(sorted(other_ab_chains)),
            'recorded_antigen_chains': ';'.join(antigen_chains),
            'false_antibody_as_antigen_chains': ';'.join(sorted(false_assignment_chains)),
            'has_false_assignment': len(false_assignment_chains) > 0,
        })

    df = pd.DataFrame(rows)
    out_path = AUDITS_DIR / "stage2_multiinstance_antigen_masking.csv"
    df.to_csv(out_path, index=False)

    logger.info(f"Wrote {out_path}")
    logger.info(f"Number of multi-INSTANCE PDBs: {len(multi_instance_pdbs)}")
    logger.info(f"Number of INSTANCEs in multi-INSTANCE PDBs (checked): {len(df)}")
    logger.info(f"Number of candidate antigen assignments involving another antibody chain: {n_false_assignments}")

    if n_false_assignments == 0:
        logger.info("RESULT: 0 false antigen assignments — PASS")
    else:
        logger.error(f"RESULT: {n_false_assignments} false antigen assignments remain — FAIL")

    return n_false_assignments == 0, len(multi_instance_pdbs), len(df), n_false_assignments


if __name__ == "__main__":
    ok, n_pdbs, n_instances, n_false = main()
    sys.exit(0 if ok else 1)
