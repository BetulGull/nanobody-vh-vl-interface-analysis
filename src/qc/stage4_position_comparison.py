"""
Stage 4 — Position-Level Amino-Acid Distributions and VHH-vs-VH Divergence
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.spatial.distance import jensenshannon
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

CANONICAL_20 = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Nonstandard/unmapped residues (e.g. 'X' for pyroglutamate/PCA -- a genuine
# post-translational N-terminal cyclization, not a data error) are tracked
# as an explicit extra category so frequencies still sum to 1.0. They
# contribute zero weight to all Kyte-Doolittle-based statistics (get_kd
# returns None for 'X'), which is the scientifically correct treatment --
# they are just not silently dropped from the frequency table itself.
NONSTANDARD_CATEGORY = 'X'
ALL_AA_CATEGORIES = CANONICAL_20 + [NONSTANDARD_CATEGORY]


def shannon_entropy(freqs: np.ndarray) -> float:
    freqs = freqs[freqs > 0]
    if len(freqs) == 0:
        return 0.0
    return float(-(freqs * np.log2(freqs)).sum())


def position_aa_distribution(long_df: pd.DataFrame, group: str) -> dict:
    """Returns {(imgt_position, insertion_code): {aa: freq, ...}} for present residues only."""
    sub = long_df[(long_df['group'] == group) & (long_df['residue_present_in_sequence'] == True)]
    result = {}
    for (pos, ins), grp in sub.groupby(['imgt_position', 'insertion_code']):
        counts = grp['amino_acid'].value_counts()
        n = counts.sum()
        freqs = {aa: counts.get(aa, 0) / n for aa in ALL_AA_CATEGORIES}
        result[(pos, ins)] = {'freqs': freqs, 'n': n, 'counts': counts.to_dict()}
    return result


def compute_position_stats(long_df: pd.DataFrame, group: str) -> pd.DataFrame:
    sub = long_df[(long_df['group'] == group) & (long_df['residue_present_in_sequence'] == True)]
    rows = []
    for (pos, ins), grp in sub.groupby(['imgt_position', 'insertion_code']):
        n = len(grp)
        counts = grp['amino_acid'].value_counts()
        freqs = np.array([counts.get(aa, 0) / n for aa in CANONICAL_20])
        rows.append({
            'imgt_position': pos, 'insertion_code': ins,
            f'n_{group}': n,
            f'entropy_{group}': shannon_entropy(freqs),
            f'mean_KD_{group}': grp['kyte_doolittle'].mean(),
            f'acidic_fraction_{group}': grp['residue_acidic'].mean(),
            f'basic_fraction_{group}': grp['residue_basic'].mean(),
            f'aromatic_fraction_{group}': grp['residue_aromatic'].mean(),
            f'positive_KD_fraction_{group}': grp['residue_positive_KD'].mean(),
        })
    return pd.DataFrame(rows)


def compute_jsd(long_df: pd.DataFrame) -> pd.DataFrame:
    vhh_dist = position_aa_distribution(long_df, 'VHH')
    vh_dist = position_aa_distribution(long_df, 'VH')

    all_positions = set(vhh_dist) | set(vh_dist)
    rows = []
    for (pos, ins) in all_positions:
        vhh_freqs = np.array([vhh_dist.get((pos, ins), {'freqs': {}})['freqs'].get(aa, 0) for aa in CANONICAL_20]) \
            if (pos, ins) in vhh_dist else np.zeros(20)
        vh_freqs = np.array([vh_dist.get((pos, ins), {'freqs': {}})['freqs'].get(aa, 0) for aa in CANONICAL_20]) \
            if (pos, ins) in vh_dist else np.zeros(20)

        if vhh_freqs.sum() == 0 or vh_freqs.sum() == 0:
            jsd = None
        else:
            jsd = float(jensenshannon(vhh_freqs, vh_freqs, base=2) ** 2)  # JS divergence (squared distance = divergence)

        rows.append({'imgt_position': pos, 'insertion_code': ins, 'JSD': jsd})

    return pd.DataFrame(rows)


def main():
    long_df = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long.parquet")

    vhh_stats = compute_position_stats(long_df, 'VHH')
    vh_stats = compute_position_stats(long_df, 'VH')
    jsd_df = compute_jsd(long_df)

    # Base position metadata (region, primary_binary_interface, tetrad, Stage-3 metrics) -- constant across sequences
    meta = long_df.drop_duplicates(subset=['imgt_position', 'insertion_code'])[
        ['imgt_position', 'insertion_code', 'region', 'primary_binary_interface', 'vhh_tetrad_position',
         'interface_frequency_5', 'mean_delta_sasa', 'median_delta_sasa', 'observation_fraction']
    ]

    comparison = meta.merge(vhh_stats, on=['imgt_position', 'insertion_code'], how='left')
    comparison = comparison.merge(vh_stats, on=['imgt_position', 'insertion_code'], how='left')
    comparison = comparison.merge(jsd_df, on=['imgt_position', 'insertion_code'], how='left')

    comparison['delta_mean_KD'] = comparison['mean_KD_VHH'] - comparison['mean_KD_VH']
    comparison['delta_acidic_fraction'] = comparison['acidic_fraction_VHH'] - comparison['acidic_fraction_VH']
    comparison['delta_basic_fraction'] = comparison['basic_fraction_VHH'] - comparison['basic_fraction_VH']
    comparison['delta_aromatic_fraction'] = comparison['aromatic_fraction_VHH'] - comparison['aromatic_fraction_VH']
    comparison['delta_positive_KD_fraction'] = comparison['positive_KD_fraction_VHH'] - comparison['positive_KD_fraction_VH']

    comparison = comparison.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)

    out_csv = PROCESSED_DIR / "stage4_position_chemistry_comparison.csv"
    out_parquet = PROCESSED_DIR / "stage4_position_chemistry_comparison.parquet"
    comparison.to_csv(out_csv, index=False)
    comparison.to_parquet(out_parquet, index=False)
    logger.info(f"Wrote {out_csv} / .parquet: {len(comparison)} positions")

    # --- AA frequency long table ---
    freq_rows = []
    for group in ['VHH', 'VH']:
        dist = position_aa_distribution(long_df, group)
        for (pos, ins), d in dist.items():
            for aa in ALL_AA_CATEGORIES:
                freq_rows.append({
                    'group': group, 'imgt_position': pos, 'insertion_code': ins,
                    'amino_acid': aa, 'count': d['counts'].get(aa, 0),
                    'frequency': d['freqs'][aa], 'n_sequences_with_position': d['n'],
                })
    freq_df = pd.DataFrame(freq_rows)
    freq_out = PROCESSED_DIR / "stage4_position_aa_frequencies.csv"
    freq_df.to_csv(freq_out, index=False)
    logger.info(f"Wrote {freq_out}: {len(freq_df)} rows")

    return comparison, freq_df


if __name__ == "__main__":
    main()
