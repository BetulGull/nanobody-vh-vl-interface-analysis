"""
Figure 4B -- CDR3 secondary structure via DSSP (mkdssp v4.6.1).

Native experimental coordinates only, no minimization/relaxation/remodeling.
Writes a minimal single-chain PDB per structure (H chain, IMGT 1-128,
VL never loaded for VH -- same convention as the rest of Figure 4),
runs mkdssp, parses the 8-state DSSP code for each CDR3 residue, and
collapses to Zhang's reported G/T/E/S/C reporting subset (Results S3.4)
-- other DSSP states (H, I, B, -) are retained in source data under
'dssp_raw_state' but not part of the Zhang-matched reporting set.
"""
import sys
import os
import subprocess
import tempfile
import time
from pathlib import Path
import pandas as pd

# mkdssp (homebrew binary) segfaults when invoked with the conda/micromamba
# env's DYLD/CONDA library-path variables inherited -- a minimal env avoids
# the conflict. Verified: full env -> returncode -11 (SIGSEGV); minimal
# PATH+HOME env -> returncode 0, correct DSSP output.
_DSSP_ENV = {'PATH': '/opt/homebrew/bin:/usr/bin:/bin', 'HOME': os.environ.get('HOME', '')}
_MKDSSP_BIN = '/opt/homebrew/bin/mkdssp'

sys.path.insert(0, str(Path(__file__).parent))
from final_manuscript_structures import get_variable_domain_atoms, write_pdb

ROOT = Path("/Users/betulgul/Desktop/staj")
CDR3_RANGE = (105, 117)
ZHANG_STATES = {'G', 'T', 'E', 'S', 'C'}


def run_dssp_one(instance, group, tmpdir):
    atoms = get_variable_domain_atoms(instance, 'H')
    pdb_path = Path(tmpdir) / f"{instance.replace('/', '_')}.pdb"
    write_pdb([('A', atoms, None)], pdb_path)

    dssp_out = Path(tmpdir) / f"{instance.replace('/', '_')}.dssp"
    try:
        result = subprocess.run([_MKDSSP_BIN, '--output-format', 'dssp', str(pdb_path), str(dssp_out)],
                                 capture_output=True, text=True, timeout=30, env=_DSSP_ENV)
        if result.returncode != 0 or not dssp_out.exists():
            return [], f"mkdssp failed: {result.stderr[:200]}"
    except Exception as e:
        return [], str(e)

    rows = []
    with open(dssp_out) as f:
        lines = f.readlines()
    started = False
    for line in lines:
        if line.startswith('  #  RESIDUE'):
            started = True
            continue
        if not started or len(line) < 17:
            continue
        try:
            resnum_str = line[5:10].strip()
            if not resnum_str:
                continue
            resnum = int(resnum_str)
            ins_code = line[10].strip()
            ss = line[16].strip()
            if not ss:
                ss = '-'
        except (ValueError, IndexError):
            continue
        if CDR3_RANGE[0] <= resnum <= CDR3_RANGE[1]:
            # DSSP's blank/unassigned state ('-') is the coil state, i.e. Zhang's
            # reported "C" -- map it explicitly rather than leaving it in OTHER.
            zhang_state = 'C' if ss == '-' else (ss if ss in ZHANG_STATES else 'OTHER')
            rows.append({'INSTANCE': instance, 'group': group, 'cdr3_position': resnum, 'cdr3_ins_code': ins_code,
                         'dssp_raw_state': ss, 'dssp_zhang_state': zhang_state})
    return rows, None


def main(limit=None, out_suffix=''):
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]
    instances = [(i, 'VHH') for i in vhh_df['INSTANCE']] + [(i, 'VH') for i in vh_df['INSTANCE']]
    if limit:
        instances = instances[:limit]

    all_rows = []
    errors = []
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, (instance, group) in enumerate(instances):
            rows, err = run_dssp_one(instance, group, tmpdir)
            all_rows.extend(rows)
            if err:
                errors.append({'INSTANCE': instance, 'group': group, 'error': err})
            if (idx + 1) % 100 == 0:
                elapsed = time.time() - t0
                rate = (idx + 1) / elapsed
                print(f"{idx+1}/{len(instances)} done, {elapsed:.1f}s, {rate:.2f}/s, ETA {(len(instances)-idx-1)/rate:.0f}s", flush=True)

    df = pd.DataFrame(all_rows)
    df.to_csv(ROOT / f"data/interim/figure4B_dssp_cdr3{out_suffix}.csv", index=False)
    err_df = pd.DataFrame(errors)
    err_df.to_csv(ROOT / f"data/interim/figure4B_dssp_errors{out_suffix}.csv", index=False)
    print(f"DONE. {len(instances)} instances, {len(errors)} errors, {len(df)} residue rows, {time.time()-t0:.1f}s")
    return df


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=None)
    p.add_argument('--suffix', type=str, default='')
    args = p.parse_args()
    main(limit=args.limit, out_suffix=args.suffix)
