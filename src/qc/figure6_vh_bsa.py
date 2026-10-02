"""
Figure 6 -- conventional-VH CDR3-to-former-interface packing BSA.

Reuses stage7_bsa.compute_pairwise_bsa VERBATIM (frozen, unmodified) --
validated in this revision to exactly reproduce (diff=0.0) the frozen VHH
nonadj_BSA_pair_total values before being trusted for a fresh VH-side
computation. Only the atom input differs: VH heavy-chain-only atoms from
the native paired VH-VL structure (VL never loaded).
"""
import sys
import time
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from final_manuscript_structures import get_variable_domain_atoms
from stage7_bsa import compute_pairwise_bsa, NONADJACENT_PACKING_TARGET

ROOT = Path("/Users/betulgul/Desktop/staj")


def main(limit=None):
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]
    instances = vh_df['INSTANCE'].tolist()
    if limit:
        instances = instances[:limit]

    rows = []
    t0 = time.time()
    for idx, inst in enumerate(instances):
        try:
            atoms = get_variable_domain_atoms(inst, 'H')
            result = compute_pairwise_bsa(atoms, NONADJACENT_PACKING_TARGET)
            row = {'INSTANCE': inst, 'group': 'VH'}
            row.update(result)
            rows.append(row)
        except Exception as e:
            rows.append({'INSTANCE': inst, 'group': 'VH', 'status': 'ERROR', 'error': str(e)})
        if (idx + 1) % 200 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            print(f"{idx+1}/{len(instances)} done, {elapsed:.1f}s, {rate:.2f}/s, ETA {(len(instances)-idx-1)/rate:.0f}s", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "data/interim/figure6_vh_bsa.csv", index=False)
    print(f"DONE. {len(out)} instances, status counts:")
    print(out['status'].value_counts())
    return out


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=None)
    args = p.parse_args()
    main(limit=args.limit)
