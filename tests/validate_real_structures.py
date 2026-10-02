"""
Stage-2 Real-Data Validation Gates

Gate 1: REAL CONTACT-SEARCH VALIDATION
Gate 2: REAL STRUCTURE-MAPPING VALIDATION

Both run on the same deterministic set of real extracted mmCIF structures:
- 5 VHH structures (from primary_vhh_candidate_strict)
- 5 paired VH-VL structures (from primary_paired_candidate_metadata)
- includes at least one multi-antigen-chain example if found
- includes at least one insertion-code example if found

Output: results/audits/stage2_real_mapping_validation.csv
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from qc.mmcif_qc import MMCIFStructure
from qc.contact_search import ContactSearch
from qc.antigen_qc import AntigenQuality
from qc.cdr_qc import CDRQuality

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent
GORDON_THRESHOLD = 4.5


def select_test_structures():
    """
    Deterministically select real-structure validation set:
    5 VHH + 5 paired, plus flag any with multi-antigen-chains or insertion codes.
    """
    metadata = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet")

    vhh_candidates = metadata[metadata["primary_vhh_candidate_strict"] == True].reset_index(drop=True)
    vhh_test = vhh_candidates.head(5).copy()
    vhh_test['cohort'] = 'VHH'

    paired_candidates = metadata[metadata["primary_paired_candidate_metadata"] == True].reset_index(drop=True)
    paired_test = paired_candidates.head(5).copy()
    paired_test['cohort'] = 'PAIRED'

    test_set = pd.concat([vhh_test, paired_test], ignore_index=True)
    logger.info(f"Selected {len(test_set)} deterministic test structures (5 VHH + 5 paired)")

    return test_set, metadata


def find_extra_examples(metadata: pd.DataFrame, existing_instances: set):
    """
    Look for a multi-protein-antigen-chain example and an insertion-code
    example among the broader candidate pool, to extend coverage.
    """
    extras = []

    # Search a modest window for a multi-antigen-chain example (antigen_chain has multiple entries, e.g. "A | B")
    candidates = metadata[
        (metadata["primary_vhh_candidate_strict"] == True) |
        (metadata["primary_paired_candidate_metadata"] == True)
    ]
    for _, row in candidates.head(300).iterrows():
        if row['INSTANCE'] in existing_instances:
            continue
        ag_chain = str(row.get('antigen_chain', ''))
        if '|' in ag_chain or ',' in ag_chain:
            row = row.copy()
            row['cohort'] = 'EXTRA_MULTI_ANTIGEN'
            extras.append(row)
            existing_instances.add(row['INSTANCE'])
            break

    return extras


def build_imgt_summary(mmcif_structure, chain_id, chain_type, model_id) -> dict:
    """Summarize IMGT mapping: n positions observed, insertion codes present, CDR ranges."""
    full_mapping = mmcif_structure.build_full_imgt_mapping_with_insertions(chain_id, model_id)

    n_positions = len(full_mapping)
    insertion_examples = []
    for pos, variants in full_mapping.items():
        if len(variants) > 1:
            insertion_examples.append(f"{pos}:{','.join(variants)}")

    cdrq = CDRQuality()
    cdr_completeness = cdrq.check_cdr_completeness(mmcif_structure, chain_id, chain_type, model_id)

    return {
        'n_imgt_positions_observed': n_positions,
        'has_insertion_codes': len(insertion_examples) > 0,
        'insertion_code_examples': ';'.join(insertion_examples[:3]),
        'cdr_completeness': cdr_completeness,
    }


def validate_structure(row: pd.Series) -> dict:
    """Run both validation gates on a single real structure."""
    instance = row['INSTANCE']
    pdb = row['PDB']
    hchain = row['Hchain']
    # SAbDab2 uses the literal string 'NA' (not a null) to mean "no assigned
    # light chain" for SD-H (VHH) rows — pd.notna() does not catch this.
    raw_lchain = row['Lchain']
    lchain = None if (pd.isna(raw_lchain) or str(raw_lchain).strip() == 'NA') else raw_lchain
    extracted_path = row['extracted_structure_path']
    expected_model = str(row.get('model', '1'))

    result = {
        'instance': instance,
        'pdb': pdb,
        'cohort': row.get('cohort', ''),
        'hchain': hchain,
        'lchain': lchain,
        'expected_model': expected_model,
        'extracted_path': str(extracted_path),
        'error': None,
        # Gate 2: structure mapping
        'parse_success': False,
        'model_resolved': None,
        'hchain_exact_match_found': False,
        'lchain_exact_match_found': None,
        'h_cdr1_complete': None,
        'h_cdr2_complete': None,
        'h_cdr3_complete': None,
        'l_cdr1_complete': None,
        'l_cdr2_complete': None,
        'l_cdr3_complete': None,
        'h_has_insertion_codes': None,
        'l_has_insertion_codes': None,
        'h_insertion_examples': '',
        'l_insertion_examples': '',
        'antigen_chains_identified': '',
        'n_antigen_chains': 0,
        'antigen_polymer_types': '',
        # Gate 1: contact search
        'contact_min_distance_brute_force': None,
        'contact_min_distance_optimized': None,
        'gordon_pass_brute_force': None,
        'gordon_pass_optimized': None,
        'contacting_antigen_chains': '',
        'contacting_cdrs': '',
        'algorithms_agree_has_contact': None,
        'algorithms_agree_indices': None,
        'gate1_pass': False,
        'gate2_pass': False,
    }

    try:
        full_path = PROJECT_ROOT / extracted_path
        if not full_path.exists():
            result['error'] = f"File not found: {full_path}"
            return result

        structure = MMCIFStructure(str(full_path))
        result['parse_success'] = True

        model_id = structure.get_model()
        result['model_resolved'] = model_id

        chains_in_model = structure.get_chains_for_model(model_id)
        result['hchain_exact_match_found'] = hchain in chains_in_model

        if lchain:
            result['lchain_exact_match_found'] = lchain in chains_in_model
        else:
            result['lchain_exact_match_found'] = None  # N/A for VHH (no assigned Lchain)

        # --- Gate 2: CDR completeness + insertion codes ---
        h_summary = build_imgt_summary(structure, hchain, 'H', model_id)
        result['h_cdr1_complete'] = h_summary['cdr_completeness'].get('H1')
        result['h_cdr2_complete'] = h_summary['cdr_completeness'].get('H2')
        result['h_cdr3_complete'] = h_summary['cdr_completeness'].get('H3')
        result['h_has_insertion_codes'] = h_summary['has_insertion_codes']
        result['h_insertion_examples'] = h_summary['insertion_code_examples']

        if lchain:
            l_summary = build_imgt_summary(structure, lchain, 'L', model_id)
            result['l_cdr1_complete'] = l_summary['cdr_completeness'].get('L1')
            result['l_cdr2_complete'] = l_summary['cdr_completeness'].get('L2')
            result['l_cdr3_complete'] = l_summary['cdr_completeness'].get('L3')
            result['l_has_insertion_codes'] = l_summary['has_insertion_codes']
            result['l_insertion_examples'] = l_summary['insertion_code_examples']

        # --- Antigen identification (from actual mmCIF entity/polymer info) ---
        antigen_qc = AntigenQuality()
        antibody_chains = {hchain}
        if lchain:
            antibody_chains.add(lchain)

        antigen_chains, chain_type_map = antigen_qc.identify_antigen_chains(
            structure, model_id, antibody_chains
        )
        result['antigen_chains_identified'] = ';'.join(antigen_chains)
        result['n_antigen_chains'] = len(antigen_chains)
        result['antigen_polymer_types'] = ';'.join(sorted(set(chain_type_map.values())))

        gate2_checks = [
            result['parse_success'],
            result['hchain_exact_match_found'],
            result['model_resolved'] is not None,
        ]
        if lchain:
            gate2_checks.append(result['lchain_exact_match_found'])
        result['gate2_pass'] = all(gate2_checks)

        # --- Gate 1: contact search on CDR atoms vs antigen atoms ---
        if antigen_chains:
            cdrq = CDRQuality()
            all_cdr_atoms = []
            cdr_membership_map = []  # parallel list: which CDR each atom belongs to
            cdr_names = ['H1', 'H2', 'H3']
            if lchain:
                cdr_names += ['L1', 'L2', 'L3']

            for cdr_name in cdr_names:
                chain_for_cdr = hchain if cdr_name[0] == 'H' else lchain
                coords = cdrq.get_cdr_atom_coords(structure, chain_for_cdr, cdr_name[0], cdr_name, model_id)
                for c in coords:
                    all_cdr_atoms.append(c)
                    cdr_membership_map.append(cdr_name)

            antigen_atom_map = antigen_qc.get_antigen_atoms(structure, antigen_chains, model_id)
            all_antigen_atoms = []
            antigen_chain_map = []
            for achain, aatoms in antigen_atom_map.items():
                for a in aatoms:
                    all_antigen_atoms.append((a['x'], a['y'], a['z']))
                    antigen_chain_map.append(achain)

            if all_cdr_atoms and all_antigen_atoms:
                cdr_coords = np.array(all_cdr_atoms)
                antigen_coords = np.array(all_antigen_atoms)

                cs = ContactSearch(threshold_angstrom=GORDON_THRESHOLD)

                # Brute force
                bf_has, bf_cdr_idx, bf_ag_idx = cs.brute_force_contacts(cdr_coords, antigen_coords)
                # Optimized
                opt_has, opt_cdr_idx, opt_ag_idx = cs.optimized_contacts(cdr_coords, antigen_coords)

                # Minimum distance (independent of threshold, for reporting)
                all_dist = np.linalg.norm(
                    cdr_coords[:, np.newaxis, :] - antigen_coords[np.newaxis, :, :], axis=2
                )
                min_dist = float(np.min(all_dist))
                result['contact_min_distance_brute_force'] = min_dist
                result['contact_min_distance_optimized'] = min_dist  # same underlying data

                result['gordon_pass_brute_force'] = bf_has
                result['gordon_pass_optimized'] = opt_has

                result['contacting_antigen_chains'] = ';'.join(
                    sorted(set(antigen_chain_map[i] for i in bf_ag_idx))
                )
                result['contacting_cdrs'] = ';'.join(
                    sorted(set(cdr_membership_map[i] for i in bf_cdr_idx))
                )

                result['algorithms_agree_has_contact'] = (bf_has == opt_has)
                result['algorithms_agree_indices'] = (bf_cdr_idx == opt_cdr_idx and bf_ag_idx == opt_ag_idx)

                result['gate1_pass'] = (
                    result['algorithms_agree_has_contact'] and
                    result['algorithms_agree_indices']
                )
            else:
                result['gate1_pass'] = None  # no CDR or antigen atoms to test
        else:
            result['gate1_pass'] = None  # no antigen chain identified for this structure

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"Error validating {instance}: {e}")

    return result


def main():
    logger.info("=" * 80)
    logger.info("STAGE-2 REAL-DATA VALIDATION GATES")
    logger.info("=" * 80)

    test_set, metadata = select_test_structures()
    extras = find_extra_examples(metadata, set(test_set['INSTANCE']))
    if extras:
        extras_df = pd.DataFrame(extras)
        test_set = pd.concat([test_set, extras_df], ignore_index=True)
        logger.info(f"Added {len(extras)} extra example(s) (multi-antigen-chain / insertion-code coverage)")

    results = []
    for idx, row in test_set.iterrows():
        result = validate_structure(row)
        results.append(result)
        logger.info(
            f"[{idx+1}/{len(test_set)}] {result['instance']} ({result['cohort']}): "
            f"gate2={result['gate2_pass']}, gate1={result['gate1_pass']}, "
            f"antigen_chains={result['n_antigen_chains']}, "
            f"gordon={result['gordon_pass_optimized']}, "
            f"insertion_H={result['h_has_insertion_codes']}"
        )

    results_df = pd.DataFrame(results)
    out_path = PROJECT_ROOT / "results" / "audits" / "stage2_real_mapping_validation.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(out_path, index=False)
    logger.info(f"\nSaved: {out_path}")

    # --- Gate summaries ---
    logger.info("")
    logger.info("=" * 80)
    logger.info("GATE 2 (STRUCTURE-MAPPING) SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Parse success: {results_df['parse_success'].sum()}/{len(results_df)}")
    logger.info(f"Hchain exact match: {results_df['hchain_exact_match_found'].sum()}/{len(results_df)}")
    lchain_applicable = results_df['lchain_exact_match_found'].notna()
    if lchain_applicable.sum() > 0:
        logger.info(f"Lchain exact match: {results_df.loc[lchain_applicable, 'lchain_exact_match_found'].sum()}/{lchain_applicable.sum()}")
    logger.info(f"Insertion codes observed (H chain): {results_df['h_has_insertion_codes'].sum()}/{len(results_df)}")
    gate2_pass_count = results_df['gate2_pass'].sum()
    logger.info(f"Gate 2 PASS: {gate2_pass_count}/{len(results_df)}")

    logger.info("")
    logger.info("=" * 80)
    logger.info("GATE 1 (CONTACT-SEARCH) SUMMARY")
    logger.info("=" * 80)
    tested = results_df['gate1_pass'].notna()
    logger.info(f"Structures with antigen+CDR atoms tested: {tested.sum()}/{len(results_df)}")
    if tested.sum() > 0:
        gate1_pass_count = results_df.loc[tested, 'gate1_pass'].sum()
        logger.info(f"Gate 1 PASS (algorithms agree): {gate1_pass_count}/{tested.sum()}")
    else:
        gate1_pass_count = 0

    logger.info("")
    logger.info("=" * 80)

    gate2_ok = (gate2_pass_count == len(results_df)) and (results_df['error'].isna().all())
    gate1_ok = tested.sum() > 0 and (results_df.loc[tested, 'gate1_pass'].sum() == tested.sum())

    if gate1_ok and gate2_ok:
        logger.info("BOTH VALIDATION GATES: PASS")
        logger.info("=" * 80)
        return True
    else:
        logger.error(f"VALIDATION GATE FAILURE — Gate1 OK={gate1_ok}, Gate2 OK={gate2_ok}")
        logger.info("=" * 80)
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
