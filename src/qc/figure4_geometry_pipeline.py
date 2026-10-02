"""
Figure 4A/4B -- shared full-cohort CDR3 geometry computation.

For every representative structure (VHH and conventional VH), extracts
IMGT-numbered heavy-chain atoms once (via the existing, validated
final_manuscript_structures / stage3_variable_domain machinery -- no new
numbering logic), then computes:

  - minimum heavy-atom distance between every real CDR3 residue and every
    residue in each of four target sets: the 15-position
    SHIELDABLE_NON_CDR3_INTERFACE, CDR1, FR2, CDR2;
  - FR2_CDR3_Dis (centroid-based, exact frozen Stage-6 definition, reused
    unmodified, applied fresh to VH coordinates).

VL is never loaded for the VH group (get_variable_domain_atoms(..., 'H')
only ever touches the H chain) -- "omitted from the calculation", not an
experimentally isolated domain.

Output: one long-format parquet per instance-batch, concatenated at the
end. Designed to run once and be reused by both Figure 4A and Figure 4B.
"""
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from final_manuscript_structures import get_variable_domain_atoms, get_row

ROOT = Path("/Users/betulgul/Desktop/staj")

SHIELDABLE_15 = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}
CDR3_RANGE = (105, 117)
CDR1_RANGE = (27, 38)
FR2_RANGE = (39, 55)
CDR2_RANGE = (56, 65)
TETRAD = {42, 49, 50, 52}

TARGET_SETS = {
    'SHIELDABLE_15': ('set', SHIELDABLE_15),
    'CDR1': ('range', CDR1_RANGE),
    'FR2': ('range', FR2_RANGE),
    'CDR2': ('range', CDR2_RANGE),
}


def group_atoms_by_position(atoms):
    """dict (imgt_pos, ins_code) -> list of atom dicts, in first-seen order (build_residue_scheme order)."""
    d = {}
    for a in atoms:
        key = (a['auth_seq_id'], a['ins_code'] if a['ins_code'] else '')
        d.setdefault(key, []).append(a)
    return d


def ordered_residue_keys(pos_dict, lo, hi):
    """Residue keys within [lo,hi], in the order build_residue_scheme naturally
    produced them (already correct IMGT insertion-order per this project's
    established convention -- same order used throughout Stage 3/4/6/7/8)."""
    seen = []
    seen_set = set()
    for a_key in pos_dict.keys():
        pos = a_key[0]
        if lo <= pos <= hi and a_key not in seen_set:
            seen.append(a_key)
            seen_set.add(a_key)
    return seen


def min_dist(atoms1, atoms2):
    p1 = np.array([[a['x'], a['y'], a['z']] for a in atoms1])
    p2 = np.array([[a['x'], a['y'], a['z']] for a in atoms2])
    diff = p1[:, None, :] - p2[None, :, :]
    d = np.sqrt((diff ** 2).sum(axis=-1))
    return float(d.min())


def centroid(atoms):
    p = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    return p.mean(axis=0)


def process_instance(instance, group):
    atoms = get_variable_domain_atoms(instance, 'H')
    pos_dict = group_atoms_by_position(atoms)

    # native, true-sequence-order CDR3 residue keys (as returned by build_residue_scheme)
    cdr3_keys_native_order = []
    seen = set()
    for a in atoms:
        key = (a['auth_seq_id'], a['ins_code'] if a['ins_code'] else '')
        if CDR3_RANGE[0] <= a['auth_seq_id'] <= CDR3_RANGE[1] and key not in seen:
            cdr3_keys_native_order.append(key)
            seen.add(key)

    n_cdr3 = len(cdr3_keys_native_order)
    internal_flags = {k: False for k in cdr3_keys_native_order}
    if n_cdr3 > 4:
        for k in cdr3_keys_native_order[2:-2]:
            internal_flags[k] = True

    # FR2_CDR3_Dis: centroid(IMGT42) to centroid of each internal CDR3 residue, min
    fr2_cdr3_dis = None
    pos42_key = (42, '')
    if pos42_key in pos_dict and n_cdr3 > 4:
        c42 = centroid(pos_dict[pos42_key])
        dists = []
        for k in cdr3_keys_native_order:
            if internal_flags[k] and k in pos_dict:
                dists.append(float(np.linalg.norm(c42 - centroid(pos_dict[k]))))
        if dists:
            fr2_cdr3_dis = min(dists)

    rows = []
    for cdr3_key in cdr3_keys_native_order:
        if cdr3_key not in pos_dict:
            continue
        cdr3_pos, cdr3_ins = cdr3_key
        cdr3_atoms = pos_dict[cdr3_key]
        for set_name, (kind, spec) in TARGET_SETS.items():
            if kind == 'set':
                target_keys = [k for k in pos_dict.keys() if k[0] in spec]
            else:
                lo, hi = spec
                target_keys = ordered_residue_keys(pos_dict, lo, hi)
            for tkey in target_keys:
                tpos, tins = tkey
                d = min_dist(cdr3_atoms, pos_dict[tkey])
                rows.append({
                    'INSTANCE': instance, 'group': group,
                    'cdr3_position': cdr3_pos, 'cdr3_ins_code': cdr3_ins,
                    'target_set': set_name, 'target_position': tpos, 'target_ins_code': tins,
                    'min_heavy_atom_distance': d,
                    'is_internal_cdr3': internal_flags[cdr3_key],
                    'tetrad': tpos in TETRAD,
                })

    return rows, fr2_cdr3_dis, n_cdr3


def main(limit=None, out_suffix=''):
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]

    instances = [(i, 'VHH') for i in vhh_df['INSTANCE']] + [(i, 'VH') for i in vh_df['INSTANCE']]
    if limit:
        instances = instances[:limit]

    all_rows = []
    fr2_rows = []
    n_fail = 0
    t0 = time.time()
    for idx, (instance, group) in enumerate(instances):
        try:
            rows, fr2_cdr3_dis, n_cdr3 = process_instance(instance, group)
            all_rows.extend(rows)
            fr2_rows.append({'INSTANCE': instance, 'group': group, 'FR2_CDR3_Dis': fr2_cdr3_dis, 'n_cdr3_residues': n_cdr3})
        except Exception as e:
            n_fail += 1
            fr2_rows.append({'INSTANCE': instance, 'group': group, 'FR2_CDR3_Dis': None, 'n_cdr3_residues': None, 'error': str(e)})
        if (idx + 1) % 50 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            print(f"{idx+1}/{len(instances)} done, {elapsed:.1f}s elapsed, {rate:.2f}/s, ETA {(len(instances)-idx-1)/rate:.0f}s", flush=True)

    dist_df = pd.DataFrame(all_rows)
    fr2_df = pd.DataFrame(fr2_rows)
    dist_df.to_parquet(ROOT / f"data/interim/figure4_pairwise_distances{out_suffix}.parquet", index=False)
    fr2_df.to_csv(ROOT / f"data/interim/figure4_fr2_cdr3_dis{out_suffix}.csv", index=False)
    print(f"DONE. {len(instances)} instances, {n_fail} failures, {len(dist_df)} distance rows, {time.time()-t0:.1f}s total")
    return dist_df, fr2_df


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=None)
    p.add_argument('--suffix', type=str, default='')
    args = p.parse_args()
    main(limit=args.limit, out_suffix=args.suffix)
