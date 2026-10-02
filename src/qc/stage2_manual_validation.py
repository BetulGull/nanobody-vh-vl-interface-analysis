"""
Stage-2 Manual Validation

Deterministically re-derives QC values for a sample of structures (>=10
strict-primary VHH, >=10 paired, plus every rare failure category with
<=10-ish examples) directly from mmCIF, independent re-computation used as a
sanity cross-check against the master table's recorded values.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.cdr_qc import CDRQuality
from qc.antigen_qc import AntigenQuality
from qc.contact_search import ContactSearch
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def recompute_and_compare(row: pd.Series, pdb_antibody_chains: dict) -> dict:
    """Independently recompute key QC values and compare to stored master-table values."""
    instance = row['INSTANCE']
    hchain = row['Hchain']
    lchain = normalize_lchain(row['Lchain'])
    extracted_path = PROJECT_ROOT / row['extracted_structure_path']

    result = {
        'instance': instance,
        'pdb': row['PDB'],
        'hchain': hchain,
        'lchain': lchain,
        'recorded_h_cdr_complete': row.get('qc_h_all_cdr_complete'),
        'recomputed_h_cdr_complete': None,
        'recorded_gordon_contact': row.get('qc_gordon_cdr_antigen_contact'),
        'recomputed_gordon_contact': None,
        'recorded_antigen_chains': row.get('antigen_chains'),
        'recomputed_antigen_chains': None,
        'match_h_cdr': None,
        'match_gordon': None,
        'match_antigen': None,
        'error': None,
    }

    try:
        s = MMCIFStructure(str(extracted_path))
        model_id = s.get_model()

        cdrq = CDRQuality()
        h_complete = cdrq.check_cdr_completeness(s, hchain, 'H', model_id)
        recomputed_h_all = all([h_complete.get('H1', False), h_complete.get('H2', False), h_complete.get('H3', False)])
        result['recomputed_h_cdr_complete'] = recomputed_h_all

        aq = AntigenQuality()
        antibody_chains = {hchain}
        if lchain:
            antibody_chains.add(lchain)
        antibody_chains |= pdb_antibody_chains.get(row['PDB'], set())
        antigen_chains, _ = aq.identify_antigen_chains(s, model_id, antibody_chains)
        result['recomputed_antigen_chains'] = ';'.join(antigen_chains)

        if not antigen_chains:
            # No antigen -> Gordon contact is unambiguously False, matching
            # the master table's default (not "unknown"/None).
            result['recomputed_gordon_contact'] = False

        if antigen_chains:
            cdr_names = ['H1', 'H2', 'H3'] + (['L1', 'L2', 'L3'] if lchain else [])
            all_cdr_atoms = []
            for cdr_name in cdr_names:
                chain_for_cdr = hchain if cdr_name[0] == 'H' else lchain
                coords = cdrq.get_cdr_atom_coords(s, chain_for_cdr, cdr_name[0], cdr_name, model_id)
                all_cdr_atoms.extend(coords)

            antigen_atom_map = aq.get_antigen_atoms(s, antigen_chains, model_id)
            all_antigen_atoms = [(a['x'], a['y'], a['z']) for atoms in antigen_atom_map.values() for a in atoms]

            if all_cdr_atoms and all_antigen_atoms:
                cs = ContactSearch(threshold_angstrom=4.5)
                has_contact, _, _ = cs.optimized_contacts(np.array(all_cdr_atoms), np.array(all_antigen_atoms))
                result['recomputed_gordon_contact'] = has_contact

        result['match_h_cdr'] = (result['recorded_h_cdr_complete'] == result['recomputed_h_cdr_complete'])
        result['match_gordon'] = (result['recorded_gordon_contact'] == result['recomputed_gordon_contact'])
        result['match_antigen'] = (
            set(str(result['recorded_antigen_chains'] or '').split(';')) ==
            set(str(result['recomputed_antigen_chains'] or '').split(';'))
        )

    except Exception as e:
        result['error'] = str(e)

    return result


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")

    pdb_antibody_chains = {}
    for _, r in master.iterrows():
        pdb = r['PDB']
        h = r['Hchain']
        l = normalize_lchain(r['Lchain'])
        pdb_antibody_chains.setdefault(pdb, set()).add(h)
        if l:
            pdb_antibody_chains[pdb].add(l)

    samples = []

    strict_primary_vhh = master[
        (master["primary_vhh_candidate_strict"] == True) &
        (master["qc_stage2_core_structure_pass"] == True) &
        (master["qc_gordon_cdr_antigen_contact"] == True)
    ]
    samples.append(strict_primary_vhh.sample(min(12, len(strict_primary_vhh)), random_state=42))

    paired_primary = master[
        (master["cohort"] == "PAIRED") &
        (master["qc_stage2_core_structure_pass"] == True) &
        (master["qc_gordon_cdr_antigen_contact"] == True)
    ]
    samples.append(paired_primary.sample(min(12, len(paired_primary)), random_state=42))

    failure_df = pd.read_csv(PROJECT_ROOT / "results" / "tables" / "stage2_failure_reasons.csv")
    edge_case_rows = []
    for reason in failure_df["failure_reason"].unique():
        instances = failure_df[failure_df["failure_reason"] == reason]["INSTANCE"].unique()
        all_inspected = len(instances) <= 12
        n_sample = len(instances) if all_inspected else min(15, len(instances))
        sampled_instances = pd.Series(instances).sample(n_sample, random_state=42) if not all_inspected else instances
        subset = master[master["INSTANCE"].isin(sampled_instances)]
        samples.append(subset)
        edge_case_rows.append({
            "failure_reason": reason, "n_examples": len(instances),
            "n_sampled": n_sample, "all_inspected": all_inspected,
        })

    sample_df = pd.concat(samples, ignore_index=True).drop_duplicates(subset="INSTANCE")
    logger.info(f"Manual validation sample: {len(sample_df)} unique structures")

    results = []
    for _, row in sample_df.iterrows():
        r = recompute_and_compare(row, pdb_antibody_chains)
        results.append(r)

    results_df = pd.DataFrame(results)
    results_df.to_csv(AUDITS_DIR / "stage2_manual_validation.csv", index=False)

    edge_cases_df = pd.DataFrame(edge_case_rows)
    edge_cases_df.to_csv(AUDITS_DIR / "stage2_edge_cases.csv", index=False)

    n_errors = results_df["error"].notna().sum()
    n_h_mismatch = (results_df["match_h_cdr"] == False).sum()
    n_gordon_mismatch = (results_df["match_gordon"] == False).sum()
    n_antigen_mismatch = (results_df["match_antigen"] == False).sum()

    logger.info(f"\nManual validation results:")
    logger.info(f"  Total sampled: {len(results_df)}")
    logger.info(f"  Errors during recomputation: {n_errors}")
    logger.info(f"  H-CDR completeness mismatches: {n_h_mismatch}")
    logger.info(f"  Gordon contact mismatches: {n_gordon_mismatch}")
    logger.info(f"  Antigen chain mismatches: {n_antigen_mismatch}")

    all_match = (n_errors == 0 and n_h_mismatch == 0 and n_gordon_mismatch == 0 and n_antigen_mismatch == 0)
    logger.info(f"\nMANUAL VALIDATION: {'PASS' if all_match else 'FAIL — systematic bug likely'}")

    return all_match


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
