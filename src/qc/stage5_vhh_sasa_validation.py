"""
Stage 5 — VHH Native Extraction + SASA Validation Gate (spec §10)

Validates on a deterministic, diverse real-structure set before any full
Stage-5 run. Reuses the SAME validated extraction (stage3_variable_domain)
and SASA (stage3_sasa) infrastructure as Stage 3 -- no new engine.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa
from qc.antigen_qc import AntigenQuality

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

CDR3_START, CDR3_END = 105, 117


def validate_case(instance: str, master: pd.DataFrame, category: str) -> dict:
    row = master[master['INSTANCE'] == instance].iloc[0]
    hchain = row['Hchain']
    path = PROJECT_ROOT / row['extracted_structure_path']

    s = MMCIFStructure(str(path))
    model_id = s.get_model()

    # Total atoms in original Hchain (for contamination comparison)
    total_h_atoms = sum(
        1 for r in s.data.get('_atom_site', [])
        if r.get('_atom_site.auth_asym_id') == hchain and
           r.get('_atom_site.pdbx_PDB_model_num') == model_id and
           r.get('_atom_site.type_symbol', '').strip().upper() not in ('H', 'D')
    )

    vhh_atoms = extract_variable_domain_atoms(s, hchain, model_id)

    # Antigen exclusion check
    aq = AntigenQuality()
    antigen_chains, _ = aq.identify_antigen_chains(s, model_id, {hchain})
    antigen_chains_in_extraction = set(a['auth_asym_id'] for a in vhh_atoms) & set(antigen_chains)

    # CDR3 intact check: all real observed residues 105-117 present in extraction
    scheme = s.build_residue_scheme(hchain, model_id)
    cdr3_real_observed = [r for r in scheme if r['observed'] and r['auth_seq_id'] is not None
                           and CDR3_START <= r['auth_seq_id'] <= CDR3_END]
    cdr3_extracted_keys = set((a['auth_seq_id'], a['ins_code']) for a in vhh_atoms
                               if CDR3_START <= a['auth_seq_id'] <= CDR3_END)
    cdr3_scheme_keys = set((r['auth_seq_id'], r['ins_code']) for r in cdr3_real_observed)
    cdr3_intact = cdr3_scheme_keys.issubset(cdr3_extracted_keys)

    # Constant-domain exclusion check
    max_extracted_pos = max((a['auth_seq_id'] for a in vhh_atoms), default=0)
    total_real_residues = len(scheme)
    has_constant_domain = total_real_residues > 128

    # Coordinate identity / reproducibility: re-extract, compare
    vhh_atoms_2 = extract_variable_domain_atoms(s, hchain, model_id)
    coords_identical = (
        len(vhh_atoms) == len(vhh_atoms_2) and
        all((a['x'], a['y'], a['z']) == (b['x'], b['y'], b['z']) for a, b in zip(vhh_atoms, vhh_atoms_2))
    )

    # FreeSASA reproducibility
    per_res_1, total_1 = calc_sasa(vhh_atoms)
    per_res_2, total_2 = calc_sasa(vhh_atoms)
    sasa_reproducible = abs(total_1 - total_2) < 1e-9

    # Insertion codes preserved
    insertion_positions = sorted(set(f"{a['auth_seq_id']}{a['ins_code']}" for a in vhh_atoms if a['ins_code']))

    return {
        'category': category, 'INSTANCE': instance, 'Hchain': hchain,
        'model_id': model_id,
        'total_atoms_original_Hchain': total_h_atoms,
        'atoms_retained_vhh_variable_domain': len(vhh_atoms),
        'antigen_atoms_absent': len(antigen_chains_in_extraction) == 0,
        'antigen_chains_identified': ';'.join(antigen_chains),
        'has_constant_domain_in_full_chain': has_constant_domain,
        'constant_domain_excluded': max_extracted_pos <= 128,
        'cdr3_real_residues_observed': len(cdr3_real_observed),
        'cdr3_fully_intact_in_extraction': cdr3_intact,
        'insertion_coded_positions_retained': ';'.join(insertion_positions),
        'coords_identical_across_calls': coords_identical,
        'total_sasa': total_1,
        'sasa_reproducible': sasa_reproducible,
    }


def main():
    vhh_master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")

    cases = [
        ('pdb_00009imx-C', 'lowest_coverage_short_or_unresolved'),
        ('pdb_00003k3q-A', 'longest_coverage_long_CDR3_insertion'),
        ('pdb_00006zbp-FFF', 'unusual_multichar_chain_id'),
        ('pdb_00006ul6-B1', 'unusual_multichar_chain_id_2'),
        ('pdb_00008qf4-E', 'highest_resolution_1.02A'),
        ('pdb_00001zmy-A', 'lowest_resolution_3.0A'),
        ('pdb_00001jtp-A', 'multiple_antigen_chain_complex'),
        ('pdb_00007m1h-G', 'low_coverage_2'),
        ('pdb_00003zkx-C', 'low_coverage_3'),
        ('pdb_00006oca-D', 'long_CDR3_2'),
        ('pdb_00009gv3-E', 'long_CDR3_3'),
    ]

    results = [validate_case(inst, vhh_master, cat) for inst, cat in cases]
    df = pd.DataFrame(results)

    out_path = AUDITS_DIR / "stage5_vhh_sasa_validation.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}: {len(df)} rows")
    logger.info(df[['category', 'INSTANCE', 'antigen_atoms_absent', 'constant_domain_excluded',
                     'cdr3_fully_intact_in_extraction', 'coords_identical_across_calls',
                     'sasa_reproducible']].to_string(index=False))

    gate_pass = (
        df['antigen_atoms_absent'].all() and
        df['constant_domain_excluded'].all() and
        df['cdr3_fully_intact_in_extraction'].all() and
        df['coords_identical_across_calls'].all() and
        df['sasa_reproducible'].all()
    )
    logger.info(f"\nVHH SASA VALIDATION GATE: {'PASS' if gate_pass else 'FAIL'}")
    return gate_pass, df


if __name__ == "__main__":
    ok, _ = main()
    sys.exit(0 if ok else 1)
