"""
Final figures -- build the combined model-input table (Section 5-6 support).

Merges the frozen Stage-11 canonical model table (stage11_model_input.parquet,
1962 rows) with the regional sequence charge/hydropathy extension
(final_regional_sequence_physicochemistry.parquet) built in
final_regional_physicochemistry.py. No new upstream science; a join only.

Note: stage11_model_input's own `cdr3_sequence_length` column is populated for
VHH only (NaN for all VH rows) -- see final_figure_data_lineage.csv. The
regional-physicochemistry table's `CDR3_length` (aggregated fresh from the
frozen Stage-4 per-position table) covers both groups and is the column used
for the CDR3-length x group comparative extension.
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")


def main():
    s11 = pd.read_parquet(ROOT / "data/processed/stage11_model_input.parquet")
    reg = pd.read_parquet(ROOT / "data/processed/final_regional_sequence_physicochemistry.parquet")

    key_cols = {'INSTANCE', 'group', 'sequence_group_id', 'SABDAB_ID'}
    reg_cols = [c for c in reg.columns if c not in key_cols]

    merged = s11.merge(reg[['INSTANCE'] + reg_cols], on='INSTANCE', how='inner')
    assert len(merged) == len(s11) == 1962, f"unexpected row count after merge: {len(merged)}"
    assert merged['CDR3_length'].isna().sum() == 0

    merged.to_parquet(ROOT / "data/processed/final_figures_model_input.parquet", index=False)
    merged.to_csv(ROOT / "data/processed/final_figures_model_input.csv", index=False)
    print(f"Wrote {len(merged)} rows, {merged.shape[1]} columns")
    return merged


if __name__ == '__main__':
    main()
