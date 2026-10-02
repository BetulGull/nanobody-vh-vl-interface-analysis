"""
Real-Structure Validation for Stage 2

Tests contact-search and mmCIF parsing on actual extracted mmCIF files before full cohort run.
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

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent


def select_test_structures():
    """
    Deterministically select 10 test structures from our cohorts.
    5 VHH + 5 paired.
    """
    metadata = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet")

    # Select first 5 strict-primary VHH
    vhh_candidates = metadata[metadata["primary_vhh_candidate_strict"] == True].reset_index(drop=True)
    vhh_test = vhh_candidates.head(5)

    # Select first 5 paired primary
    paired_candidates = metadata[metadata["primary_paired_candidate_metadata"] == True].reset_index(drop=True)
    paired_test = paired_candidates.head(5)

    test_set = pd.concat([vhh_test, paired_test], ignore_index=True)

    logger.info(f"Selected {len(test_set)} test structures")
    logger.info(f"  VHH: {len(vhh_test)}")
    logger.info(f"  Paired: {len(paired_test)}")

    return test_set


def validate_structure(row: pd.Series) -> dict:
    """
    Validate a single structure.

    Returns:
        Dict with validation results
    """
    instance = row['INSTANCE']
    pdb = row['PDB']
    hchain = row['Hchain']
    lchain = row['Lchain']
    extracted_path = row['extracted_structure_path']

    result = {
        'instance': instance,
        'pdb': pdb,
        'hchain': hchain,
        'lchain': lchain,
        'extracted_path': str(extracted_path),
        'parse_success': False,
        'hchain_found': False,
        'lchain_found': False if pd.notna(lchain) else None,
        'has_antigen': False,
        'contact_search_tested': False,
        'contact_search_brute_force_ok': False,
        'contact_search_optimized_ok': False,
        'algorithms_agree': False,
        'error': None,
    }

    try:
        # Parse mmCIF
        full_path = PROJECT_ROOT / extracted_path
        if not full_path.exists():
            result['error'] = f"File not found: {full_path}"
            return result

        structure = MMCIFStructure(str(full_path))
        result['parse_success'] = True

        # Get model
        model_id = structure.get_model()

        # Check chains
        chains_in_model = structure.get_chains_for_model(model_id)
        result['hchain_found'] = hchain in chains_in_model

        if pd.notna(lchain):
            result['lchain_found'] = lchain in chains_in_model

        # Test antigen identification
        antigen_qc = AntigenQuality()
        antibody_chains = {hchain}
        if pd.notna(lchain):
            antibody_chains.add(lchain)

        antigen_chains, _ = antigen_qc.identify_antigen_chains(
            structure, model_id, antibody_chains
        )
        result['has_antigen'] = len(antigen_chains) > 0

        # Test contact search on actual structure
        # Get H-chain atoms
        hchain_atoms = []
        for row_atom in structure.data.get('_atom_site', []):
            if (row_atom.get('_atom_site.auth_asym_id') == hchain and
                row_atom.get('_atom_site.pdbx_PDB_model_num') == model_id):
                try:
                    hchain_atoms.append({
                        'x': float(row_atom.get('_atom_site.Cartn_x', 0)),
                        'y': float(row_atom.get('_atom_site.Cartn_y', 0)),
                        'z': float(row_atom.get('_atom_site.Cartn_z', 0)),
                    })
                except (ValueError, TypeError):
                    pass

        # Get antigen atoms
        antigen_atoms = []
        for achain in antigen_chains:
            for row_atom in structure.data.get('_atom_site', []):
                if (row_atom.get('_atom_site.auth_asym_id') == achain and
                    row_atom.get('_atom_site.pdbx_PDB_model_num') == model_id):
                    try:
                        antigen_atoms.append({
                            'x': float(row_atom.get('_atom_site.Cartn_x', 0)),
                            'y': float(row_atom.get('_atom_site.Cartn_y', 0)),
                            'z': float(row_atom.get('_atom_site.Cartn_z', 0)),
                        })
                    except (ValueError, TypeError):
                        pass

        if hchain_atoms and antigen_atoms:
            result['contact_search_tested'] = True

            # Convert to numpy arrays
            hchain_coords = np.array([[a['x'], a['y'], a['z']] for a in hchain_atoms])
            antigen_coords = np.array([[a['x'], a['y'], a['z']] for a in antigen_atoms])

            # Test both algorithms
            cs = ContactSearch(threshold_angstrom=4.5)
            bf_result = cs.brute_force_contacts(hchain_coords, antigen_coords)
            opt_result = cs.optimized_contacts(hchain_coords, antigen_coords)

            result['contact_search_brute_force_ok'] = bf_result[0] is not None
            result['contact_search_optimized_ok'] = opt_result[0] is not None

            # Check agreement
            result['algorithms_agree'] = (
                bf_result[0] == opt_result[0] and
                bf_result[1] == opt_result[1] and
                bf_result[2] == opt_result[2]
            )

    except Exception as e:
        result['error'] = str(e)
        logger.error(f"Error validating {instance}: {e}")

    return result


def main():
    """Run validation on test set."""
    logger.info("=" * 80)
    logger.info("REAL-STRUCTURE VALIDATION FOR STAGE 2")
    logger.info("=" * 80)
    logger.info("")

    test_set = select_test_structures()

    results = []
    for idx, row in test_set.iterrows():
        result = validate_structure(row)
        results.append(result)
        logger.info(f"[{idx+1}/{len(test_set)}] {result['instance']}: parse_ok={result['parse_success']}, "
                   f"hchain_found={result['hchain_found']}, antigen_found={result['has_antigen']}, "
                   f"contact_ok={result['contact_search_tested']}, agree={result['algorithms_agree']}")

    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 80)

    results_df = pd.DataFrame(results)
    logger.info(f"Parse success: {results_df['parse_success'].sum()}/{len(results_df)}")
    logger.info(f"Hchain found: {results_df['hchain_found'].sum()}/{len(results_df)}")
    logger.info(f"Lchain found: {results_df['lchain_found'].sum()}/{results_df['lchain_found'].notna().sum()}")
    logger.info(f"Antigen found: {results_df['has_antigen'].sum()}/{len(results_df)}")
    logger.info(f"Contact search tested: {results_df['contact_search_tested'].sum()}/{len(results_df)}")
    logger.info(f"Algorithms agree: {results_df['algorithms_agree'].sum()}/{results_df['contact_search_tested'].sum()}")

    # Save results
    results_df.to_csv(PROJECT_ROOT / "results" / "audits" / "stage2_real_mapping_validation.csv", index=False)
    logger.info(f"Results saved to results/audits/stage2_real_mapping_validation.csv")

    # Fail if critical validation failed
    if results_df['parse_success'].sum() < len(results_df) * 0.9:
        logger.error("VALIDATION FAILED: <90% parse success")
        return False

    if results_df['hchain_found'].sum() < len(results_df) * 0.9:
        logger.error("VALIDATION FAILED: <90% hchain found")
        return False

    if results_df['contact_search_tested'].sum() > 0:
        if results_df['algorithms_agree'].sum() < results_df['contact_search_tested'].sum():
            logger.error("VALIDATION FAILED: Algorithms disagree on contact results")
            return False

    logger.info("VALIDATION PASSED")
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
