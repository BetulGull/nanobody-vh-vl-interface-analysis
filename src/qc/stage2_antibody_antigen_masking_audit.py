"""
Stage-2 Final Semantic Audit — Antibody-Chain Antigen Masking Semantics

For every chain currently suppressed from antigen candidacy by the
multi-instance mask (because it is H/L for another SAbDab2 INSTANCE sharing
the same PDB), determine WHY it was suppressed and whether that suppression
is actually correct, using:
  - SAbDab2's own INSTANCE-specific antigen_chain / antigen_type annotation
    (strong cross-check, not sole truth source)
  - mmCIF entity/polymer identity (already confirmed protein-polymer, since
    only such chains are ever candidates in the first place)
  - SABDAB_ID relationship to the "other" instance (same sequence == same
    antibody component, not a distinct entity)

Output: results/audits/stage2_antibody_antigen_mask_semantics.csv
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.antigen_qc import AntigenQuality

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


ANTI_IDIOTYPE_STRONG_KEYWORDS = [
    'ANTI-IDIOTYPE', 'ANTI IDIOTYPE', 'IDIOTOPE', 'ANTI-ID ', 'ANTIIDIOTYPE',
    'PROTEIN A', 'PROTEIN G', 'PROTEIN L', 'FC RECEPTOR', 'FCRN', 'NEONATAL FC RECEPTOR',
    'RHEUMATOID FACTOR',
]

# Any of these appearing in antigen_name is evidence the "antigen" is simply
# another antibody/nanobody bound to the SAME real target (a sibling in a
# multi-Ab complex) or an apo crystallographic partner — NOT a deliberate
# antibody-as-antigen relationship. SAbDab2's own antigen_chain field turned
# out to list such siblings indiscriminately (see audit narrative), so this
# keyword alone must NOT be treated as confirming legitimacy.
SIBLING_ANTIBODY_KEYWORDS = [
    'NANOBODY', 'NB1', 'NB2', 'NB0', ' FAB ', 'SCFV', 'VHH', 'SINGLE-DOMAIN',
    'SDAB', 'CDR1-3 GRAFT', 'IMMUNOGLOBULIN', 'HYBRIDOMA', 'MONOCLONAL ANTIBODY',
    ' MAB ', 'MAB)', '(LIGHT CHAIN)', '(HEAVY CHAIN)', 'ANTIBODY LA2', 'ANTIBODY FRAGMENT',
]


def classify_suppressed_chain(row, suppressed_chain, other_instance_rows, pdb_antibody_chains) -> tuple:
    """
    Returns (classification, evidence_note).
    """
    own_antigen_chain_field = str(row.get('antigen_chain', '') or '')
    own_antigen_chain_list = [c.strip() for c in own_antigen_chain_field.replace('|', ',').split(',') if c.strip()]
    antigen_name = str(row.get('antigen_name', '') or '').upper()

    named_by_own_metadata = suppressed_chain in own_antigen_chain_list

    if named_by_own_metadata:
        has_strong_antiidiotype_evidence = any(kw in antigen_name for kw in ANTI_IDIOTYPE_STRONG_KEYWORDS)
        has_sibling_antibody_signal = any(kw in antigen_name for kw in SIBLING_ANTIBODY_KEYWORDS)

        if has_strong_antiidiotype_evidence and not has_sibling_antibody_signal:
            return ('POTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN',
                    f"SAbDab2 antigen_chain names '{suppressed_chain}' AND antigen_name contains a strong "
                    f"anti-idiotype/antibody-binding-protein signal with no sibling-antibody contamination "
                    f"(antigen_type={row.get('antigen_type')}, antigen_name={row.get('antigen_name')})")
        else:
            # SAbDab2's antigen_chain field names this chain, but the evidence
            # is a known-unreliable pattern for multi-INSTANCE PDBs (see
            # narrative): apo antibody-dimer self-reference, or a sibling
            # nanobody/Fab bound to the same real target that SAbDab2 lumps
            # into the antigen_name/antigen_chain fields indiscriminately.
            # Fall through to chain-relationship evidence below instead of
            # trusting this field alone.
            reason_suffix = (
                f"[own antigen_chain field DOES name '{suppressed_chain}', but antigen_name "
                f"('{row.get('antigen_name')}') shows sibling-antibody contamination "
                f"({has_sibling_antibody_signal}) or lacks a strong anti-idiotype signal "
                f"({not has_strong_antiidiotype_evidence}) — SAbDab2's own field is unreliable here, not treated as sole truth]"
            )
    else:
        reason_suffix = ''

    # Which other INSTANCE(s) own this chain as their H or L?
    owning_rows = other_instance_rows[
        (other_instance_rows['Hchain'] == suppressed_chain) |
        (other_instance_rows['Lchain'] == suppressed_chain)
    ]

    if len(owning_rows) == 0:
        return ('UNRESOLVED', 'Suppressed chain not found as H/L of any other tracked INSTANCE in this PDB (unexpected)')

    other_row = owning_rows.iloc[0]

    # Evidence 2: same SABDAB_ID as this instance -> literally the same
    # antibody sequence, i.e. a duplicate/component of the same entity, not
    # an independent chain at all.
    if pd.notna(row.get('SABDAB_ID')) and row.get('SABDAB_ID') == other_row.get('SABDAB_ID'):
        return ('SAME_ANTIBODY_COMPONENT',
                f"Owning INSTANCE {other_row['INSTANCE']} shares SABDAB_ID {row.get('SABDAB_ID')} with this row "
                f"(same antibody sequence/entity, not a distinct chain)")

    # Evidence 3: does the OTHER instance's own antigen_chain field point
    # back at THIS instance's H/L, WITH a strong anti-idiotype signal and no
    # sibling-antibody contamination in ITS antigen_name? Mutual naming alone
    # is NOT sufficient — the pdb_00004pfe case shows two apo-dimer-adjacent
    # antibody copies mutually (and wrongly) naming each other via the exact
    # same unreliable heuristic; mutuality of an unreliable field does not
    # make it reliable.
    other_antigen_field = str(other_row.get('antigen_chain', '') or '')
    other_antigen_list = [c.strip() for c in other_antigen_field.replace('|', ',').split(',') if c.strip()]
    other_antigen_name = str(other_row.get('antigen_name', '') or '').upper()
    own_chains = {row['Hchain']}
    if normalize_lchain(row['Lchain']):
        own_chains.add(normalize_lchain(row['Lchain']))

    if own_chains & set(other_antigen_list):
        other_has_strong = any(kw in other_antigen_name for kw in ANTI_IDIOTYPE_STRONG_KEYWORDS)
        other_has_sibling_signal = any(kw in other_antigen_name for kw in SIBLING_ANTIBODY_KEYWORDS)
        if other_has_strong and not other_has_sibling_signal:
            return ('POTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN',
                    f"Owning INSTANCE {other_row['INSTANCE']}'s own antigen_chain field names this INSTANCE's "
                    f"chain(s) {own_chains & set(other_antigen_list)}, AND its antigen_name "
                    f"('{other_row.get('antigen_name')}') carries a strong anti-idiotype signal with no "
                    f"sibling-antibody contamination")
        reason_suffix += (
            f" [mutual antigen_chain cross-reference exists with {other_row['INSTANCE']} but its antigen_name "
            f"('{other_row.get('antigen_name')}') also shows sibling-antibody contamination or lacks a strong "
            f"anti-idiotype signal — same unreliable pattern, not treated as confirming]"
        )

    # Default: ordinary co-deposited antibody copy targeting its own
    # (typically different) antigen in the same crystal/PDB entry.
    return ('SIBLING_CRYSTAL_COPY_OR_PARALLEL_ANTIBODY_INSTANCE',
            f"Owning INSTANCE {other_row['INSTANCE']} (SABDAB_ID={other_row.get('SABDAB_ID')}) is an ordinary "
            f"co-deposited antibody with its own antigen_chain={other_row.get('antigen_chain')}; no reliable "
            f"evidence linking it to this INSTANCE as antigen.{reason_suffix}")


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

    pdb_instance_counts = master.groupby('PDB')['INSTANCE'].nunique()
    multi_instance_pdbs = set(pdb_instance_counts[pdb_instance_counts > 1].index)

    aq = AntigenQuality()
    results = []
    n_errors = 0

    multi_rows = master[master['PDB'].isin(multi_instance_pdbs)]
    logger.info(f"Scanning {len(multi_rows)} INSTANCEs across {len(multi_instance_pdbs)} multi-INSTANCE PDBs...")

    for _, row in multi_rows.iterrows():
        if not row['qc_structure_resolved'] or not row['qc_hchain_exists']:
            continue

        pdb = row['PDB']
        hchain = row['Hchain']
        lchain = normalize_lchain(row['Lchain'])
        own_chains = {hchain}
        if lchain:
            own_chains.add(lchain)

        path = PROJECT_ROOT / row['extracted_structure_path']
        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
        except Exception:
            n_errors += 1
            continue

        # Baseline: what would be identified as antigen using ONLY own-chain exclusion
        baseline_chains, _ = aq.identify_antigen_chains(s, model_id, own_chains)

        # Current (post-fix): full multi-instance mask applied
        current_chains = set(str(row.get('antigen_chains', '') or '').split(';')) - {''}

        suppressed = set(baseline_chains) - current_chains
        if not suppressed:
            continue

        other_instance_rows = master[(master['PDB'] == pdb) & (master['INSTANCE'] != row['INSTANCE'])]

        for chain in suppressed:
            classification, evidence = classify_suppressed_chain(row, chain, other_instance_rows, pdb_antibody_chains)
            results.append({
                'INSTANCE': row['INSTANCE'],
                'PDB': pdb,
                'own_H': hchain,
                'own_L': lchain or '',
                'suppressed_chain': chain,
                'classification': classification,
                'evidence': evidence,
                'this_instance_antigen_chain_field': row.get('antigen_chain'),
                'this_instance_antigen_type_field': row.get('antigen_type'),
            })

    df = pd.DataFrame(results)
    out_path = AUDITS_DIR / "stage2_antibody_antigen_mask_semantics.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path} ({len(df)} suppressed-chain classifications, {n_errors} parse errors)")

    if len(df) > 0:
        counts = df['classification'].value_counts()
        logger.info("\nClassification breakdown:")
        logger.info(counts.to_string())

        n_legit = (df['classification'] == 'POTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN').sum()
        logger.info(f"\nPOTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN cases: {n_legit}")
        if n_legit > 0:
            logger.info("Examples:")
            logger.info(df[df['classification'] == 'POTENTIAL_LEGITIMATE_ANTIBODY_ANTIGEN'].to_string(index=False))
    else:
        logger.info("No suppressed chains found among multi-INSTANCE PDBs.")

    return df


if __name__ == "__main__":
    main()
