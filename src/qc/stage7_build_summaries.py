"""
Stage 7 — Build per-VHH, per-position-pair, target-position, and CDR3-position
summary tables, plus the combined Stage 4-7 join, from the long interaction
table and BSA results.
"""
import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from stage7_bsa import NONADJACENT_PACKING_TARGET, FULL_SHIELDABLE_TARGET

ROOT = Path("/Users/betulgul/Desktop/staj")
TABLES_DIR = ROOT / "results/tables"


def load_bsa():
    rows = []
    with open(ROOT / "data/interim/stage7/bsa_results.jsonl") as f:
        for line in f:
            rows.append(json.loads(line))
    return pd.DataFrame(rows)


def load_status():
    rows = []
    with open(ROOT / "data/interim/stage7/cohort_status.jsonl") as f:
        for line in f:
            rows.append(json.loads(line))
    return pd.DataFrame(rows)


def build_per_vhh_summary(long_df, bsa_df, status_df):
    # primary = only include_primary_noncovalent records (nonadjacent target, specific, noncovalent)
    primary = long_df[long_df['include_primary_noncovalent']]
    # sensitivity = specific noncovalent records over FULL_SHIELDABLE_TARGET (includes 118, excludes covalent)
    sensitivity = long_df[long_df['is_specific_interaction'] & ~long_df['is_covalent']]

    def agg(df, label):
        if df.empty:
            return pd.DataFrame(columns=['INSTANCE', f'{label}_n_specific_interactions_atomlevel',
                                          f'{label}_n_unique_cdr3_target_respairs'])
        g = df.groupby('INSTANCE')
        n_atom = g.size().rename(f'{label}_n_specific_interactions_atomlevel')
        n_pairs = df.groupby('INSTANCE').apply(
            lambda d: d[['cdr3_imgt_position', 'cdr3_ins_code', 'target_imgt_position']].drop_duplicates().shape[0]
        ).rename(f'{label}_n_unique_cdr3_target_respairs')
        return pd.concat([n_atom, n_pairs], axis=1).reset_index()

    primary_agg = agg(primary, 'primary')
    sens_agg = agg(sensitivity, 'sensitivity')

    summary = bsa_df.merge(primary_agg, on='INSTANCE', how='left').merge(sens_agg, on='INSTANCE', how='left')
    summary = summary.merge(status_df[['INSTANCE', 'status', 'n_records']], on='INSTANCE', how='left')

    for col in ['primary_n_specific_interactions_atomlevel', 'primary_n_unique_cdr3_target_respairs',
                'sensitivity_n_specific_interactions_atomlevel', 'sensitivity_n_unique_cdr3_target_respairs']:
        summary[col] = summary[col].fillna(0).astype(int)

    # normalized densities (NA, not zero, when denominator undefined)
    cdr3_len = None  # will merge from stage6 geometry below in combined table; keep raw counts here
    summary['primary_interactions_per_100A2_BSA'] = summary.apply(
        lambda r: (r['primary_n_unique_cdr3_target_respairs'] / r['nonadj_BSA_pair_total'] * 100)
        if pd.notna(r['nonadj_BSA_pair_total']) and r['nonadj_BSA_pair_total'] > 0 else np.nan, axis=1)

    summary.to_csv(TABLES_DIR / "stage7_per_vhh_packing_summary.csv", index=False)
    print(f"Per-VHH summary: {len(summary)} rows")
    return summary


def build_position_pair_map(long_df):
    primary = long_df[long_df['include_primary_noncovalent']]
    pair_counts = primary.groupby(['cdr3_imgt_position', 'target_imgt_position'])['INSTANCE'].nunique().reset_index()
    pair_counts.columns = ['cdr3_imgt_position', 'target_imgt_position', 'n_instances_with_interaction']
    n_total_instances = long_df['INSTANCE'].nunique()
    pair_counts['interaction_frequency'] = pair_counts['n_instances_with_interaction'] / n_total_instances
    pair_counts = pair_counts.sort_values('interaction_frequency', ascending=False)
    pair_counts.to_csv(TABLES_DIR / "stage7_position_pair_interaction_map.csv", index=False)
    print(f"Position-pair map: {len(pair_counts)} pairs")
    return pair_counts


def build_target_and_cdr3_summaries(long_df):
    primary = long_df[long_df['include_primary_noncovalent']]
    n_total = long_df['INSTANCE'].nunique()

    target_summary = primary.groupby('target_imgt_position').agg(
        n_instances_involved=('INSTANCE', 'nunique'),
        n_atom_records=('INSTANCE', 'size'),
    ).reset_index()
    target_summary['involvement_frequency'] = target_summary['n_instances_involved'] / n_total
    target_summary.to_csv(TABLES_DIR / "stage7_target_position_summary.csv", index=False)

    cdr3_summary = primary.groupby('cdr3_imgt_position').agg(
        n_instances_involved=('INSTANCE', 'nunique'),
        n_atom_records=('INSTANCE', 'size'),
    ).reset_index()
    cdr3_summary['involvement_frequency'] = cdr3_summary['n_instances_involved'] / n_total
    cdr3_summary.to_csv(TABLES_DIR / "stage7_cdr3_position_summary.csv", index=False)

    print(f"Target-position summary: {len(target_summary)} rows; CDR3-position summary: {len(cdr3_summary)} rows")
    return target_summary, cdr3_summary


def build_combined_table(summary):
    stage6 = pd.read_parquet(ROOT / "data/processed/stage6_cdr3_shielding_geometry_combined.parquet")
    stage4 = pd.read_csv(ROOT / "data/processed/stage4_sequence_interface_chemistry.csv")
    stage4 = stage4.rename(columns={'representative_INSTANCE': 'INSTANCE'})

    combined = stage6.merge(summary, on='INSTANCE', how='left', suffixes=('', '_s7'))
    if 'PDB_s7' in combined.columns:
        combined = combined.drop(columns=['PDB_s7'])
    combined = combined.merge(stage4, on='INSTANCE', how='left', suffixes=('', '_s4'))

    combined.to_csv(TABLES_DIR / "stage7_stage4to7_combined.csv", index=False)
    print(f"Combined Stage 4-7 table: {len(combined)} rows, {len(combined.columns)} columns")
    return combined


if __name__ == '__main__':
    long_df = pd.read_csv(TABLES_DIR / "stage7_interaction_records_long.csv")
    bsa_df = load_bsa()
    status_df = load_status()
    summary = build_per_vhh_summary(long_df, bsa_df, status_df)
    build_position_pair_map(long_df)
    build_target_and_cdr3_summaries(long_df)
    build_combined_table(summary)
