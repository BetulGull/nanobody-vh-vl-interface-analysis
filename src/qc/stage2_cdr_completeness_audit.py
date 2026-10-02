"""
Stage-2 Closure Audit — CDR Completeness Semantics

Produces targeted, human-readable examples proving the distinction between:
  A. IMGT numbering placeholder gaps (no real residue — not a completeness failure)
  B. Real residues with missing coordinates (genuine completeness failure)
  C. Insertion-coded CDR residues (111, 111A, 111B kept distinct)
  D. A fully complete, unremarkable CDR

Output: results/audits/stage2_cdr_completeness_semantics.csv
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.cdr_qc import CDRQuality

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def describe_cdr(row: pd.Series, chain_type: str, cdr_name: str, category: str) -> dict:
    hchain = row['Hchain']
    lchain = None if (pd.isna(row['Lchain']) or str(row['Lchain']).strip() == 'NA') else row['Lchain']
    chain_id = hchain if chain_type == 'H' else lchain
    extracted_path = PROJECT_ROOT / row['extracted_structure_path']

    s = MMCIFStructure(str(extracted_path))
    model_id = s.get_model()
    cdrq = CDRQuality()

    scheme = s.build_residue_scheme(chain_id, model_id)
    cdr_start, cdr_end = cdrq.get_cdr_range(cdr_name)
    real_residues = cdrq._find_cdr_real_residues(scheme, cdr_start, cdr_end)

    completeness = cdrq.check_cdr_completeness(s, chain_id, chain_type, model_id)
    missing = cdrq.get_missing_cdr_positions(s, chain_id, chain_type, model_id)

    if real_residues is None:
        real_residues = []

    actual_seq = ''.join(_three_to_one(r['mon_id']) for r in real_residues)
    imgt_positions = [
        f"{r['auth_seq_id']}{r['ins_code']}" if r['observed'] else "unresolved"
        for r in real_residues
    ]
    numbering_gaps = [
        pos for pos in range(cdr_start, cdr_end + 1)
        if pos not in {r['auth_seq_id'] for r in real_residues if r['observed'] and r['auth_seq_id'] is not None}
    ]
    observed_residues = [f"{r['auth_seq_id']}{r['ins_code']}" for r in real_residues if r['observed']]
    missing_real = [f"label_seq_id={m['label_seq_id']}({m['mon_id']}, {m['reason']})" for m in missing.get(cdr_name, [])]

    return {
        'category': category,
        'PDB': row['PDB'],
        'INSTANCE': row['INSTANCE'],
        'chain': chain_id,
        'CDR': cdr_name,
        'actual_sequence_residues': actual_seq,
        'n_real_residues_in_cdr': len(real_residues),
        'imgt_positions_of_real_residues': ';'.join(imgt_positions),
        'unused_imgt_numbering_slots_in_range': ';'.join(str(p) for p in numbering_gaps),
        'observed_coordinate_residues': ';'.join(observed_residues),
        'missing_real_residues': ';'.join(missing_real) if missing_real else '(none)',
        'final_completeness_status': completeness.get(cdr_name),
    }


def _three_to_one(three):
    table = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
        'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V',
    }
    return table.get(three, 'X')


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")

    examples = []

    # A. Short CDR with legitimate IMGT numbering gaps but complete coordinates
    #    (pdb_00001bzq-K CDR-H1, confirmed earlier: 27-30 + 35-38 present, 31-34
    #    are unused IMGT slots for this loop's actual length, all observed complete)
    row_a = master[master['INSTANCE'] == 'pdb_00001bzq-K'].iloc[0]
    examples.append(describe_cdr(row_a, 'H', 'H1', 'A_legitimate_IMGT_gap_complete_coords'))

    # B. Genuine missing-coordinate CDR residue (found during closure-audit scan)
    row_b = master[master['INSTANCE'] == 'pdb_00001kxt-B'].iloc[0]
    examples.append(describe_cdr(row_b, 'H', 'H1', 'B_genuine_missing_coordinate_residue'))

    # C. Insertion-coded CDR residue (pdb_00001bzq-K CDR-H3 has an insertion at 112)
    examples.append(describe_cdr(row_a, 'H', 'H3', 'C_insertion_coded_residue'))

    # D. Fully complete, unremarkable CDR (H2 of the same well-resolved 1BZQ K chain)
    examples.append(describe_cdr(row_a, 'H', 'H2', 'D_fully_complete_normal_cdr'))

    df = pd.DataFrame(examples)
    out_path = AUDITS_DIR / "stage2_cdr_completeness_semantics.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}")
    logger.info(df.to_string(index=False))

    return df


if __name__ == "__main__":
    main()
