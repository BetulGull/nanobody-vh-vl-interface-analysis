"""
Stage 2 Orchestrator — Structure-Level QC and Final Structural Cohort Construction

Loads Stage-1 resolved metadata, processes every candidate structure through:
  structure resolution -> model resolution -> chain verification ->
  variable-domain / CDR completeness -> protein-antigen identification ->
  Gordon <=4.5A CDR-antigen contact -> (paired only) VH-VL physical-pair check

Produces per-row QC flags and failure reasons, checkpointed for resumability.

Does NOT perform: SASA/delta-SASA, CDR3 geometry, interface derivation,
Arpeggio, APBS, aggregation proxies, or any MD. Coordinate-level QC only.
"""

import logging
import traceback
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd
import numpy as np
import yaml

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.structure_resolver import StructureResolver
from qc.mmcif_qc import MMCIFStructure, mmCIFParseError
from qc.cdr_qc import CDRQuality
from qc.antigen_qc import AntigenQuality
from qc.contact_search import ContactSearch
from qc.pairing_qc import PairingQuality

PROJECT_ROOT = Path(__file__).parent.parent.parent
CONFIG_PATH = PROJECT_ROOT / "configs" / "stage2.yaml"
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage2" / "checkpoints"
CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / "stage2_structure_qc.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def load_config() -> Dict[str, Any]:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_input_metadata():
    logger.info("Loading Stage-1 resolved metadata...")
    resolved = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet")

    vhh_df = resolved[resolved["vhh_candidate_structural_qc"] == True].copy().reset_index(drop=True)
    logger.info(f"VHH structural-QC pool: {len(vhh_df)} rows")

    paired_df = resolved[resolved["primary_paired_candidate_metadata"] == True].copy().reset_index(drop=True)
    logger.info(f"Paired primary pool: {len(paired_df)} rows")

    return vhh_df, paired_df


def normalize_lchain(raw_lchain) -> Any:
    """SAbDab2 uses literal string 'NA' (not null) for 'no assigned light chain'."""
    if pd.isna(raw_lchain):
        return None
    if str(raw_lchain).strip() == 'NA':
        return None
    return raw_lchain


def build_pdb_antibody_chain_map(vhh_df: pd.DataFrame, paired_df: pd.DataFrame) -> Dict[str, set]:
    """
    Map PDB -> set of ALL chain IDs used as Hchain/Lchain by ANY SAbDab2
    INSTANCE in our input pools for that PDB.

    Required because a single deposited PDB entry commonly contains multiple
    independent antibody copies (multiple SAbDab2 INSTANCEs). Antigen
    identification for one INSTANCE must exclude every other INSTANCE's
    antibody chains too — not just its own H/L pair — or another antibody's
    variable-domain chains get misidentified as a "protein antigen" merely
    for being a non-self protein polymer chain in the same structure file.
    """
    chain_map: Dict[str, set] = {}
    for df in (vhh_df, paired_df):
        for _, row in df.iterrows():
            pdb = row['PDB']
            h = row['Hchain']
            l = normalize_lchain(row['Lchain'])
            chain_map.setdefault(pdb, set()).add(h)
            if l:
                chain_map[pdb].add(l)
    return chain_map


def process_one_structure(
    row: pd.Series,
    resolver: StructureResolver,
    cohort: str,
    pdb_antibody_chains: Dict[str, set],
) -> Dict[str, Any]:
    """
    Run full Stage-2 QC on a single candidate row.

    Returns a dict of QC flags / values to be merged into the master table.
    Never raises for row-local structural issues; only genuinely unexpected
    (systematic) errors propagate, which the caller treats specially.
    """
    instance = row['INSTANCE']
    pdb = row['PDB']
    hchain = row['Hchain']
    lchain = normalize_lchain(row['Lchain'])
    extracted_path = row['extracted_structure_path']
    archive_sha256 = row.get('archive_sha256', '')

    out = {
        'INSTANCE': instance,
        'cohort': cohort,
        'qc_structure_resolved': False,
        'qc_model_resolved': False,
        'resolved_model_id': None,
        'qc_hchain_exists': False,
        'qc_lchain_exists': False,
        'qc_no_assigned_light_chain': lchain is None,
        'qc_h_cdr1_complete': False,
        'qc_h_cdr2_complete': False,
        'qc_h_cdr3_complete': False,
        'qc_l_cdr1_complete': None if lchain is None else False,
        'qc_l_cdr2_complete': None if lchain is None else False,
        'qc_l_cdr3_complete': None if lchain is None else False,
        'qc_h_all_cdr_complete': False,
        'qc_l_all_cdr_complete': None if lchain is None else False,
        'qc_h3_geometry_anchor_complete': False,
        'qc_vh_vl_physically_paired': False,
        'vh_vl_min_distance_angstrom': None,
        'qc_has_candidate_protein_antigen': False,
        'antigen_chains': '',
        'n_antigen_chains': 0,
        'qc_gordon_cdr_antigen_contact': False,
        'contacting_cdrs': '',
        'contacting_antigen_chains': '',
        'qc_stage2_core_structure_pass': False,
        'qc_stage2_gordon_complex_pass': False,
        'failure_reasons': [],
        'case_collision_resolved': False,
        'parse_error': None,
    }

    # 1. Structure resolution
    resolution = resolver.resolve_structure_path(extracted_path, pdb, instance, archive_sha256)
    if not resolution['resolved']:
        out['failure_reasons'].append(resolution['error_reason'] or 'FAIL_STRUCTURE_RESOLUTION')
        out['failure_reasons'] = ';'.join(out['failure_reasons'])
        return out

    out['qc_structure_resolved'] = True
    out['case_collision_resolved'] = resolution['case_collision_resolved']
    structure_path = resolution['structure_path']

    # 2. mmCIF parse + model resolution
    try:
        structure = MMCIFStructure(structure_path)
    except mmCIFParseError as e:
        out['parse_error'] = str(e)
        out['failure_reasons'].append('FAIL_MMCIF_PARSE')
        out['failure_reasons'] = ';'.join(out['failure_reasons'])
        return out

    model_id = structure.get_model()
    out['qc_model_resolved'] = True
    out['resolved_model_id'] = model_id

    chains_in_model = structure.get_chains_for_model(model_id)

    # 3. Chain verification (case-sensitive)
    out['qc_hchain_exists'] = hchain in chains_in_model
    if not out['qc_hchain_exists']:
        out['failure_reasons'].append('FAIL_HCHAIN_NOT_FOUND')

    if lchain is not None:
        out['qc_lchain_exists'] = lchain in chains_in_model
        if not out['qc_lchain_exists']:
            out['failure_reasons'].append('FAIL_LCHAIN_NOT_FOUND')
    else:
        out['qc_lchain_exists'] = None  # not applicable for VHH

    if not out['qc_hchain_exists'] or (lchain is not None and not out['qc_lchain_exists']):
        out['failure_reasons'] = ';'.join(out['failure_reasons'])
        return out

    # 4. CDR completeness
    cdrq = CDRQuality()
    h_completeness = cdrq.check_cdr_completeness(structure, hchain, 'H', model_id)
    out['qc_h_cdr1_complete'] = h_completeness.get('H1', False)
    out['qc_h_cdr2_complete'] = h_completeness.get('H2', False)
    out['qc_h_cdr3_complete'] = h_completeness.get('H3', False)
    out['qc_h_all_cdr_complete'] = all([out['qc_h_cdr1_complete'], out['qc_h_cdr2_complete'], out['qc_h_cdr3_complete']])
    out['qc_h3_geometry_anchor_complete'] = cdrq.check_h3_geometry_anchor(structure, hchain, model_id)

    if not out['qc_h_all_cdr_complete']:
        out['failure_reasons'].append('FAIL_H_CDR_INCOMPLETE')

    if lchain is not None:
        l_completeness = cdrq.check_cdr_completeness(structure, lchain, 'L', model_id)
        out['qc_l_cdr1_complete'] = l_completeness.get('L1', False)
        out['qc_l_cdr2_complete'] = l_completeness.get('L2', False)
        out['qc_l_cdr3_complete'] = l_completeness.get('L3', False)
        out['qc_l_all_cdr_complete'] = all([out['qc_l_cdr1_complete'], out['qc_l_cdr2_complete'], out['qc_l_cdr3_complete']])
        if not out['qc_l_all_cdr_complete']:
            out['failure_reasons'].append('FAIL_L_CDR_INCOMPLETE')
    else:
        # No assigned light chain (VHH) — L-CDR flags are not applicable, not failures.
        out['qc_l_cdr1_complete'] = None
        out['qc_l_cdr2_complete'] = None
        out['qc_l_cdr3_complete'] = None
        out['qc_l_all_cdr_complete'] = None

    # 5. Antigen identification
    # Exclude this INSTANCE's own H/L chains AND every other SAbDab2
    # INSTANCE's antibody chains in the same PDB (see build_pdb_antibody_chain_map) —
    # a co-deposited second antibody copy must never be misidentified as antigen.
    aq = AntigenQuality()
    antibody_chains = {hchain}
    if lchain is not None:
        antibody_chains.add(lchain)
    antibody_chains |= pdb_antibody_chains.get(pdb, set())
    antigen_chains, chain_type_map = aq.identify_antigen_chains(structure, model_id, antibody_chains)
    out['qc_has_candidate_protein_antigen'] = len(antigen_chains) > 0
    out['antigen_chains'] = ';'.join(antigen_chains)
    out['n_antigen_chains'] = len(antigen_chains)

    if not out['qc_has_candidate_protein_antigen']:
        out['failure_reasons'].append('NO_PROTEIN_ANTIGEN_IDENTIFIED')

    # 6. Gordon CDR-antigen contact search (only over CDR-restricted atoms)
    if antigen_chains:
        all_cdr_atoms = []
        cdr_membership = []
        cdr_names = ['H1', 'H2', 'H3']
        if lchain is not None:
            cdr_names += ['L1', 'L2', 'L3']

        for cdr_name in cdr_names:
            chain_for_cdr = hchain if cdr_name[0] == 'H' else lchain
            coords = cdrq.get_cdr_atom_coords(structure, chain_for_cdr, cdr_name[0], cdr_name, model_id)
            for c in coords:
                all_cdr_atoms.append(c)
                cdr_membership.append(cdr_name)

        antigen_atom_map = aq.get_antigen_atoms(structure, antigen_chains, model_id)
        all_antigen_atoms = []
        antigen_chain_membership = []
        for achain, aatoms in antigen_atom_map.items():
            for a in aatoms:
                all_antigen_atoms.append((a['x'], a['y'], a['z']))
                antigen_chain_membership.append(achain)

        if all_cdr_atoms and all_antigen_atoms:
            cs = ContactSearch(threshold_angstrom=4.5)
            cdr_coords = np.array(all_cdr_atoms)
            antigen_coords = np.array(all_antigen_atoms)
            has_contact, cdr_idx, ag_idx = cs.optimized_contacts(cdr_coords, antigen_coords)

            out['qc_gordon_cdr_antigen_contact'] = has_contact
            if has_contact:
                out['contacting_cdrs'] = ';'.join(sorted(set(cdr_membership[i] for i in cdr_idx)))
                out['contacting_antigen_chains'] = ';'.join(sorted(set(antigen_chain_membership[i] for i in ag_idx)))
            else:
                out['failure_reasons'].append('NO_GORDON_CONTACT')
        else:
            out['failure_reasons'].append('NO_CDR_OR_ANTIGEN_ATOMS_FOR_CONTACT_SEARCH')

    # 7. VH-VL physical-pair sanity check (paired cohort only)
    if lchain is not None:
        pq = PairingQuality(threshold_angstrom=5.0)
        vh_atoms = pq.get_variable_domain_atoms(structure, hchain, model_id)
        vl_atoms = pq.get_variable_domain_atoms(structure, lchain, model_id)
        is_paired, min_dist = pq.check_physical_pair(vh_atoms, vl_atoms)
        out['qc_vh_vl_physically_paired'] = is_paired
        out['vh_vl_min_distance_angstrom'] = min_dist if min_dist != float('inf') else None
        if not is_paired:
            out['failure_reasons'].append('FAIL_VH_VL_NOT_PHYSICALLY_PAIRED')

    # 8. Core structure pass / Gordon complex pass
    core_checks = [out['qc_hchain_exists'], out['qc_h_all_cdr_complete']]
    if lchain is not None:
        core_checks += [out['qc_lchain_exists'], out['qc_l_all_cdr_complete'], out['qc_vh_vl_physically_paired']]
    out['qc_stage2_core_structure_pass'] = all(core_checks)

    out['qc_stage2_gordon_complex_pass'] = (
        out['qc_stage2_core_structure_pass'] and
        out['qc_has_candidate_protein_antigen'] and
        out['qc_gordon_cdr_antigen_contact']
    )

    out['failure_reasons'] = ';'.join(out['failure_reasons'])
    return out


def run_cohort(
    df: pd.DataFrame,
    cohort_name: str,
    resolver: StructureResolver,
    checkpoint_prefix: str,
    pdb_antibody_chains: Dict[str, set],
) -> pd.DataFrame:
    """Process a full cohort with checkpointing every 100 rows."""
    results = []
    n = len(df)
    systematic_error_count = 0

    for idx, row in df.iterrows():
        try:
            result = process_one_structure(row, resolver, cohort_name, pdb_antibody_chains)
        except Exception as e:
            logger.error(f"UNEXPECTED error processing {row.get('INSTANCE', '?')}: {e}")
            logger.error(traceback.format_exc())
            systematic_error_count += 1
            result = {
                'INSTANCE': row['INSTANCE'],
                'cohort': cohort_name,
                'qc_structure_resolved': False,
                'failure_reasons': f'UNEXPECTED_EXCEPTION: {e}',
            }
            # Stop if errors look systematic (parser-level), not row-local
            if systematic_error_count >= 20:
                logger.error(f"STOPPING: {systematic_error_count} unexpected exceptions — likely systematic error, not row-local.")
                raise

        results.append(result)

        if (idx + 1) % 100 == 0 or (idx + 1) == n:
            partial_df = pd.DataFrame(results)
            checkpoint_path = CHECKPOINT_DIR / f"{checkpoint_prefix}_{idx+1:05d}.parquet"
            partial_df.to_parquet(checkpoint_path, index=False)
            logger.info(f"[{cohort_name}] Checkpoint: {idx+1}/{n} processed -> {checkpoint_path.name}")

    return pd.DataFrame(results)


def stage2_main():
    logger.info("=" * 80)
    logger.info("STAGE 2 — STRUCTURE-LEVEL QC AND FINAL STRUCTURAL COHORT CONSTRUCTION")
    logger.info("=" * 80)

    vhh_df, paired_df = load_input_metadata()
    resolver = StructureResolver()
    pdb_antibody_chains = build_pdb_antibody_chain_map(vhh_df, paired_df)
    logger.info(f"Built PDB->antibody-chain map for multi-instance antigen masking: {len(pdb_antibody_chains)} PDBs")

    logger.info("")
    logger.info(f"Processing VHH cohort ({len(vhh_df)} structures)...")
    vhh_qc = run_cohort(vhh_df, 'VHH', resolver, 'vhh_checkpoint', pdb_antibody_chains)

    logger.info("")
    logger.info(f"Processing paired cohort ({len(paired_df)} structures)...")
    paired_qc = run_cohort(paired_df, 'PAIRED', resolver, 'paired_checkpoint', pdb_antibody_chains)

    return vhh_df, vhh_qc, paired_df, paired_qc


if __name__ == "__main__":
    stage2_main()
