"""
Final figures — regional sequence charge/hydropathy extension (Section 6).

Derived ONLY from the frozen Stage-4 per-position, per-representative table
(stage4_interface_projection_long_all_instances.parquet), aggregated to the
frozen final 533 VHH / 1429 VH representative cohorts. No re-numbering, no
re-parsing of raw sequences.

Charge convention: source-verified against Zhang et al. 2025 Methods
(D=-1, E=-1, R=+1, K=+1, H=+0.1, pH 7.4) -- identical to the convention
already used in Stage 8's PSH/PPC/PNC machinery.
Hydropathy: mean Kyte-Doolittle over observed region residues (same KD
table used throughout this project since Stage 4).
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")

CHARGE_MAP = {'D': -1.0, 'E': -1.0, 'R': 1.0, 'K': 1.0, 'H': 0.1}

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}


def main():
    long_df = pd.read_parquet(ROOT / "data/processed/stage4_interface_projection_long_all_instances.parquet")
    long_df = long_df.rename(columns={'representative_INSTANCE': 'INSTANCE'})

    vhh_final = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")['INSTANCE']
    vh_final_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_final = vh_final_df[vh_final_df['is_primary_representative'] == True]['INSTANCE']
    final_instances = set(vhh_final) | set(vh_final)

    df = long_df[long_df['INSTANCE'].isin(final_instances)].copy()
    df = df[df['residue_present_in_sequence'] == True]  # PRESENT_OBSERVED or PRESENT_UNRESOLVED, never coordinate-only
    df['charge'] = df['amino_acid'].map(CHARGE_MAP).fillna(0.0)
    df['is_former_interface'] = df['imgt_position'].isin(FULL_EMPIRICAL_INTERFACE)
    df['region_ext'] = df.apply(lambda r: 'FULL_EMPIRICAL_INTERFACE_28' if r['is_former_interface'] else r['region'], axis=1)

    rows = []
    region_defs = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    for instance, g in df.groupby('INSTANCE'):
        group = g['group'].iloc[0]
        seq_group_id = g['sequence_group_id'].iloc[0]
        pdb = g['SABDAB_ID'].iloc[0]
        row = {'INSTANCE': instance, 'group': group, 'sequence_group_id': seq_group_id, 'SABDAB_ID': pdb}
        for region in region_defs:
            sub = g[g['region'] == region]
            row[f'{region}_charge'] = sub['charge'].sum()
            row[f'{region}_mean_hydropathy'] = sub['kyte_doolittle'].mean()
            row[f'{region}_n_residues'] = len(sub)

        cdr_sub = g[g['region'].isin(['CDR1', 'CDR2', 'CDR3'])]
        row['all_CDR_charge'] = cdr_sub['charge'].sum()
        row['all_CDR_mean_hydropathy'] = cdr_sub['kyte_doolittle'].mean()
        row['all_CDR_n_residues'] = len(cdr_sub)

        row['full_domain_charge'] = g['charge'].sum()
        row['full_domain_mean_hydropathy'] = g['kyte_doolittle'].mean()
        row['full_domain_n_residues'] = len(g)

        interface_sub = g[g['is_former_interface']]
        row['FULL_EMPIRICAL_INTERFACE_28_charge'] = interface_sub['charge'].sum()
        row['FULL_EMPIRICAL_INTERFACE_28_mean_hydropathy'] = interface_sub['kyte_doolittle'].mean()
        row['FULL_EMPIRICAL_INTERFACE_28_n_residues'] = len(interface_sub)

        row['CDR3_length'] = len(g[g['region'] == 'CDR3'])

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_parquet(ROOT / "data/processed/final_regional_sequence_physicochemistry.parquet", index=False)
    out.to_csv(ROOT / "data/processed/final_regional_sequence_physicochemistry.csv", index=False)
    print(f"Wrote {len(out)} rows ({(out.group=='VHH').sum()} VHH, {(out.group=='VH').sum()} VH)")
    print(out.groupby('group')['CDR3_length'].describe())
    return out


if __name__ == '__main__':
    main()
