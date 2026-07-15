# -*- coding: utf-8 -*-
"""
Master driver for the FOLDED-FI MBO campaign (self-contained in new_simulation/).

Phase A : evaluate the remaining initial LHS designs 12..39 using their EXISTING
          cutout geometries (build -> solve -> folded post). Designs 0..11 were
          seeded from retained ODBs (seed_existing.py).
Phase B : run the qLogEHVI BO loop (run_mbo.py) from design 40 to 400.

Resumable: re-run any time; anything already valid in results.csv is skipped.

  python new_simulation/run_campaign.py --q 3 --cpus 8   # this machine (24 cores)
  python new_simulation/run_campaign.py --q 1 --cpus 8   # q=1 machine (Abaqus 2017)
"""
import os
import sys
import shutil
import argparse
import subprocess

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from run_mbo import read_results_rows, OBJ_FI_COL
from abaqus_runner import run_batch_iteration, run_full_iteration

BASE_DIR = THIS_DIR
INIT_BUDGET = 40


def done_specs():
    _, rows = read_results_rows()
    done = set()
    for r in rows:
        try:
            if r.get(OBJ_FI_COL, '').strip() and r.get('MAX_SE', '').strip():
                done.add(int(r['SPECIMEN_NUM']))
        except (KeyError, ValueError):
            pass
    return done


def clean_sim(n):
    d = os.path.join(BASE_DIR, 'sim_%d' % n)
    if os.path.isdir(d):
        shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, default=3)
    ap.add_argument('--cpus', type=int, default=8)
    ap.add_argument('--max-retries', type=int, default=1)
    args = ap.parse_args()

    print('=' * 60)
    print('  FOLDED-FI MBO campaign  (q=%d, cpus/job=%d)' % (args.q, args.cpus))
    print('  BASE_DIR :', BASE_DIR)
    print('=' * 60)

    # ---------- Phase A: initial LHS (existing geometries in specimen/) ----------
    # Runs whatever of 0..39 is not yet in results.csv. On this machine 0..11
    # were seeded from retained ODBs, so only 12..39 run; on a fresh machine
    # (e.g. the q=1 Abaqus-2017 box) all 40 are evaluated.
    done = done_specs()
    todo = [n for n in range(0, INIT_BUDGET) if n not in done]
    print('\n[Phase A] initial LHS remaining (0..39): %s' % todo)
    for s in range(0, len(todo), args.q):
        batch = todo[s:s + args.q]
        print('\n--- initial batch %s ---' % batch)
        for n in batch:
            clean_sim(n)
        statuses = run_batch_iteration(batch, BASE_DIR, cpus_per_job=args.cpus)
        for n in batch:
            ok = statuses[n][0]
            if ok:
                continue
            print('  RETRY initial %d' % n)
            clean_sim(n)
            for _ in range(args.max_retries):
                ok2, _p, _r = run_full_iteration(n, BASE_DIR, cpus=args.cpus)
                if ok2:
                    break
                clean_sim(n)

    # verify initial 40 complete
    done = done_specs()
    missing = [n for n in range(0, INIT_BUDGET) if n not in done]
    if missing:
        print('\nWARNING: initial designs still missing: %s' % missing)
        print('Re-run run_campaign.py to retry them before BO relies on them.')
    else:
        print('\n[Phase A] all 40 initial designs evaluated.')

    # ---------- Phase B: BO loop 40..399 ----------
    print('\n[Phase B] launching qLogEHVI BO loop (run_mbo.py)...')
    rc = subprocess.call([sys.executable, os.path.join(BASE_DIR, 'run_mbo.py'),
                          '--q', str(args.q), '--cpus', str(args.cpus)])
    print('\ncampaign finished (BO rc=%d).' % rc)


if __name__ == '__main__':
    main()
