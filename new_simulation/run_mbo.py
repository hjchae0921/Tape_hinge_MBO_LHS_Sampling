# -*- coding: utf-8 -*-
"""
Main MBO driver for the tape-hinge Bezier cutout problem.

Spec (from MBO.md, with the q-batch update; budget extended to 270 in a
later run by raising `TOTAL_BUDGET` from 240 to 270):
  - 270 total designs = 40 initial LHS (already in results.csv) + 230 BO
  - Objectives (maximize):  y0 = 1 - MAX_FI_ALL,  y1 = MAX_SE
  - Reference point     :  (0.0, 400.0)
  - Kernel              :  Matern12 (nu=0.5), fixed noise variance 1e-6
  - Hyperparameters     :  MAP refit every iteration
                           (lengthscale ~ Gamma(3, 6), outputscale ~ Gamma(2, 0.15)
                            -- BoTorch SingleTaskGP default priors)
  - Acquisition         :  qLogEHVI (BoTorch recommended log-domain
                           numerical fix of qLogEHVI; arXiv:2310.20708.
                           Same acquisition mathematically, more stable
                           optimization. Reported value is log-domain.)
  - Batch / parallelism :  q = 3, each solver job uses 4 cpus
                           (build/post are sequential to avoid contention
                            on Abaqus CAE licenses and on results.csv)
  - Random seed         :  42  (deterministic for reproduction)
  - Abaqus 2024 Explicit

Run:
    python mbo/run_mbo.py                       # auto-loop until total budget
    python mbo/run_mbo.py --q 3 --cpus 4        # batch of 3, 4 cpus per solver
    python mbo/run_mbo.py --max-iters 6         # smoke test (cap BO count this run)
    python mbo/run_mbo.py --total-budget 270    # override total design count

Resume:
    Just re-run.  The driver reads results.csv, derives next SPECIMEN_NUM
    from `max(SPECIMEN_NUM) + 1`, and continues until `--total-budget` is
    reached.  Raising `--total-budget` between runs simply extends the BO
    tail (same seed, same GP refit, same q=3 qLogEHVI loop).

Failure handling:
    A specimen that fails the batch step is retried individually (single
    sequential pass). If the retry also fails, a placeholder row is
    appended to results.csv (only SPECIMEN_NUM filled) so the loop does
    not stall, and `mbo_log.csv` records status='failed'.
"""
import os
import sys
import csv
import math
import time
import shutil
import argparse
import datetime

import numpy as np
import torch

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
# This campaign is SELF-CONTAINED in new_simulation/: results.csv, specimen/,
# sim_<N>/, plot/ and the Abaqus scripts all live alongside this driver.
BASE_DIR = THIS_DIR
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from bezier_cutout import DV_NAMES, DV_BOUNDS, write_cutout_csv
from gp_model import (
    build_gp, propose_next, hypervolume, hyperparams_summary, REF_POINT,
)
from abaqus_runner import run_batch_iteration, run_full_iteration
from plotter import plot_pareto


# -----------------------------------------------------------------
# Constants
# -----------------------------------------------------------------
SEED           = 42
TOTAL_BUDGET   = 400                          # MBO.md original: 240 (40 LHS + 200 BO).
                                              # Extended to 270 to add 30 more BO specimens
                                              # on top of the first 200.
INIT_BUDGET    = 40
BO_BUDGET      = TOTAL_BUDGET - INIT_BUDGET   # 230
BATCH_Q        = 3
CPUS_PER_JOB   = 4

SPECIMEN_DIR   = os.path.join(BASE_DIR, 'specimen')
RESULTS_CSV    = os.path.join(BASE_DIR, 'results.csv')
MBO_LOG_CSV    = os.path.join(BASE_DIR, 'mbo_log.csv')
PLOT_DIR       = os.path.join(BASE_DIR, 'plot')

# results.csv stores `Eccentricity`, but DV_NAMES uses `ecc`.
# These keys are in the canonical DV order (matches DV_NAMES).
DV_KEYS_IN_CSV = ['y1', 'x1', 'Dx', 'y3', 'yP2', 'fP3', 'Eccentricity']

# Objective FI is now the FULLY-FOLDED FI (MAX_FI_FOLDED). MAX_FI_ALL etc.
# are kept for reference. 17-column canonical header.
RESULTS_HEADER_16 = [
    'SPECIMEN_NUM',
    'MAX_FI_ALL', 'MAX_FI_PINCHING', 'MAX_FI_FOLDING', 'MAX_FI_FOLDED', 'MAX_SE',
    'Eccentricity', 'y1', 'x1', 'Dx', 'y3', 'yP2', 'fP3',
    'x2_derived', 'x3_derived', 'xP3_derived', 'yP3_derived',
]
OBJ_FI_COL = 'MAX_FI_FOLDED'      # objective: y0 = 1 - fully-folded FI


# -----------------------------------------------------------------
# results.csv helpers
# -----------------------------------------------------------------
def _trim_to_canonical_header(header, rows):
    """Drop trailing empty columns left over from manual PC merges."""
    n_keep = len(RESULTS_HEADER_16)
    # If header already matches, no-op.
    if header[:n_keep] == RESULTS_HEADER_16:
        return RESULTS_HEADER_16, [r[:n_keep] for r in rows]
    return header, rows


def read_results_rows():
    """Return (header, [dict rows]) from results.csv (sorted by SPECIMEN_NUM)."""
    if not os.path.isfile(RESULTS_CSV):
        return RESULTS_HEADER_16, []
    with open(RESULTS_CSV, 'r', newline='') as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        return RESULTS_HEADER_16, []
    hdr = rows[0]
    data = rows[1:]
    hdr, data = _trim_to_canonical_header(hdr, data)
    dict_rows = []
    for r in data:
        if not r or not r[0].strip().isdigit():
            continue
        # pad short rows
        r = r + [''] * (len(hdr) - len(r))
        dict_rows.append(dict(zip(hdr, r)))
    dict_rows.sort(key=lambda d: int(d['SPECIMEN_NUM']))
    return hdr, dict_rows


def append_placeholder_row(specimen_num):
    """Append a row with only SPECIMEN_NUM populated, so the index advances."""
    hdr, rows = read_results_rows()
    placeholder = {c: '' for c in RESULTS_HEADER_16}
    placeholder['SPECIMEN_NUM'] = str(specimen_num)
    rows.append(placeholder)
    rows.sort(key=lambda d: int(d['SPECIMEN_NUM']))
    with open(RESULTS_CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=RESULTS_HEADER_16,
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def rows_to_tensors(rows):
    """Filter rows with valid (1-FI, SE); return (X, Y) double tensors."""
    X_list, Y_list = [], []
    for r in rows:
        try:
            x_vals = []
            for k in DV_KEYS_IN_CSV:
                s = r[k].strip()
                if not s:
                    raise ValueError('missing DV %s' % k)
                x_vals.append(float(s))
            fi_s = r[OBJ_FI_COL].strip()      # fully-folded FI objective
            se_s = r['MAX_SE'].strip()
            if not fi_s or not se_s:
                continue
            fi = float(fi_s)
            se = float(se_s)
            if math.isnan(fi) or math.isnan(se):
                continue
        except (KeyError, ValueError):
            continue
        X_list.append(x_vals)
        Y_list.append([1.0 - fi, se])
    if not X_list:
        return (torch.zeros((0, len(DV_NAMES)), dtype=torch.double),
                torch.zeros((0, 2), dtype=torch.double))
    X = torch.tensor(X_list, dtype=torch.double)
    Y = torch.tensor(Y_list, dtype=torch.double)
    return X, Y


def get_bounds_tensor():
    return torch.tensor(DV_BOUNDS.T, dtype=torch.double)   # (2, 7)


def next_specimen_index(rows):
    """Resume-aware: 40 if init-only, otherwise max+1."""
    if not rows:
        return INIT_BUDGET
    return max(int(r['SPECIMEN_NUM']) for r in rows) + 1


# -----------------------------------------------------------------
# mbo_log.csv helpers
# -----------------------------------------------------------------
def append_mbo_log(specimen_num, status, dv, hv, acq_val, hp_str, note=''):
    write_header = not os.path.isfile(MBO_LOG_CSV)
    with open(MBO_LOG_CSV, 'a', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        if write_header:
            w.writerow(['timestamp', 'SPECIMEN_NUM', 'status',
                        'HV', 'qLogEHVI']
                       + DV_NAMES
                       + ['gp_hyperparams', 'note'])
        ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        dv_round = [('%.6f' % v) for v in dv]
        w.writerow([ts, specimen_num, status,
                    '%.6f' % hv, '%.6e' % acq_val]
                   + dv_round
                   + [hp_str, note])


# -----------------------------------------------------------------
# Main
# -----------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description='Tape-hinge MBO driver')
    ap.add_argument('--q', type=int, default=BATCH_Q,
                    help='qLogEHVI batch size (default: %d)' % BATCH_Q)
    ap.add_argument('--cpus', type=int, default=CPUS_PER_JOB,
                    help='CPUs per Abaqus solver job (default: %d)'
                         % CPUS_PER_JOB)
    ap.add_argument('--total-budget', type=int, default=TOTAL_BUDGET,
                    help='Total design count cap (initial LHS + BO). '
                         'Driver stops once next SPECIMEN_NUM reaches this. '
                         '(default: %d)' % TOTAL_BUDGET)
    ap.add_argument('--max-iters', type=int, default=BO_BUDGET,
                    help='Cap BO specimen count THIS run (default: %d). '
                         'Use to slice a long extension into shorter sessions; '
                         '--total-budget controls the overall cap.'
                         % BO_BUDGET)
    ap.add_argument('--max-retries', type=int, default=1,
                    help='Retries per failing specimen (default: 1)')
    return ap.parse_args()


def _check_no_collision(indices):
    """Abort if any sim_<N>/ or specimen/cutout_points_<N>.csv already exists."""
    for n in indices:
        sim_dir = os.path.join(BASE_DIR, 'sim_%d' % n)
        csv_path = os.path.join(SPECIMEN_DIR, 'cutout_points_%d.csv' % n)
        if os.path.isdir(sim_dir):
            print('ABORT: sim_%d already exists. Resolve manually before resuming.'
                  % n)
            sys.exit(2)
        if os.path.isfile(csv_path):
            print('ABORT: %s already exists. Resolve manually before resuming.'
                  % csv_path)
            sys.exit(2)


def main():
    args = parse_args()

    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.use_deterministic_algorithms(False)

    bounds_t = get_bounds_tensor()

    total_budget = args.total_budget
    bo_budget    = total_budget - INIT_BUDGET

    print('=' * 64)
    print('  Tape-hinge MBO  (qLogEHVI, Matern12, noise_var=1e-6, seed=42)')
    print('  Reference point :', REF_POINT)
    print('  Total budget    : %d  (init %d + BO %d)'
          % (total_budget, INIT_BUDGET, bo_budget))
    print('  Batch size q    : %d   cpus/solver-job : %d'
          % (args.q, args.cpus))
    print('  results.csv     :', RESULTS_CSV)
    print('  mbo_log.csv     :', MBO_LOG_CSV)
    print('=' * 64)

    bo_done = 0
    iter_idx = 0
    while bo_done < args.max_iters:
        iter_idx += 1
        # ---- reload state (resume-safe) -----------------------------
        _, rows = read_results_rows()
        X, Y = rows_to_tensors(rows)
        next_idx = next_specimen_index(rows)
        if next_idx >= total_budget:
            print('Total budget reached (next idx %d >= %d).'
                  % (next_idx, total_budget))
            break

        q_eff = min(args.q,
                    total_budget - next_idx,
                    args.max_iters - bo_done)
        if q_eff <= 0:
            break

        cand_indices = list(range(next_idx, next_idx + q_eff))
        _check_no_collision(cand_indices)

        print('\n' + '=' * 64)
        print('  BO iter %d  q_eff=%d  SPECIMEN_NUM=%s  valid data n=%d'
              % (iter_idx, q_eff, cand_indices, X.shape[0]))
        print('=' * 64)

        if X.shape[0] == 0:
            print('ABORT: no valid rows in results.csv -- cannot fit GP.')
            sys.exit(2)

        # ---- fit GP, propose q candidates jointly -------------------
        t0 = time.time()
        model = build_gp(X, Y, bounds_t)
        cand_X, acq_val = propose_next(
            model, Y, bounds_t, seed=SEED + iter_idx, q=q_eff,
        )
        hp_str = hyperparams_summary(model)
        dv_batch = cand_X.cpu().numpy()             # (q_eff, 7)
        hv_before = hypervolume(Y)
        t_bo = time.time() - t0

        print('  BO step in %.1fs   joint qLogEHVI=%.4e   HV(before)=%.4f'
              % (t_bo, acq_val, hv_before))
        for i, n in enumerate(cand_indices):
            print('  cand #%d  (SPECIMEN_NUM=%d):' % (i + 1, n))
            for n_, v_ in zip(DV_NAMES, dv_batch[i]):
                print('    %-5s = %.6f' % (n_, v_))

        # ---- write cutout CSVs --------------------------------------
        for n, dv in zip(cand_indices, dv_batch):
            csv_path = os.path.join(SPECIMEN_DIR,
                                    'cutout_points_%d.csv' % n)
            write_cutout_csv(csv_path, dv)
        print('  cutout CSVs written for %s' % cand_indices)

        # ---- run batch: build seq, solver parallel, post seq --------
        statuses = run_batch_iteration(
            cand_indices, BASE_DIR, cpus_per_job=args.cpus,
        )

        # ---- per-specimen retry for failures (sequential) -----------
        for n in cand_indices:
            ok, phase, rcs = statuses[n]
            if ok:
                continue
            print('  RETRY specimen %d (batch failed at phase=%s)' % (n, phase))
            sim_dir = os.path.join(BASE_DIR, 'sim_%d' % n)
            if os.path.isdir(sim_dir):
                shutil.rmtree(sim_dir, ignore_errors=True)
            retry_ok = False
            for attempt in range(args.max_retries):
                ok2, phase2, rcs2 = run_full_iteration(
                    n, BASE_DIR, cpus=args.cpus,
                )
                if ok2:
                    statuses[n] = (True, 'done', rcs2)
                    retry_ok = True
                    break
                print('    retry %d still failed at phase=%s rcs=%s'
                      % (attempt + 1, phase2, rcs2))
                if os.path.isdir(sim_dir):
                    shutil.rmtree(sim_dir, ignore_errors=True)
            if not retry_ok:
                # keep last failure info in statuses[n]
                statuses[n] = (False, phase2 if args.max_retries else phase,
                               rcs2 if args.max_retries else rcs)

        # ---- log results + handle placeholders ----------------------
        _, rows_after = read_results_rows()
        Y_after_full = rows_to_tensors(rows_after)[1]
        hv_after = hypervolume(Y_after_full)

        for i, n in enumerate(cand_indices):
            ok, phase, rcs = statuses[n]
            dv_i = dv_batch[i]
            if not ok:
                append_placeholder_row(n)
                append_mbo_log(n, 'failed', dv_i,
                               hv_before, acq_val, hp_str,
                               note='last_phase=%s rcs=%s' % (phase, rcs))
                continue
            match = [r for r in rows_after
                     if r.get('SPECIMEN_NUM', '').strip() == str(n)]
            valid = False
            if match:
                try:
                    fi = float(match[0][OBJ_FI_COL])
                    se = float(match[0]['MAX_SE'])
                    valid = not (math.isnan(fi) or math.isnan(se))
                except (KeyError, ValueError):
                    valid = False
            if valid:
                append_mbo_log(n, 'ok', dv_i, hv_after, acq_val, hp_str)
            else:
                append_mbo_log(n, 'invalid_result', dv_i,
                               hv_before, acq_val, hp_str,
                               note='row present but MAX_FI_ALL/MAX_SE missing')

        # ---- Pareto scatter plot ------------------------------------
        try:
            Y_init_list, Y_bo_list = [], []
            for r in rows_after:
                try:
                    n_ = int(r['SPECIMEN_NUM'])
                    fi_s = r[OBJ_FI_COL].strip()
                    se_s = r['MAX_SE'].strip()
                    if not fi_s or not se_s:
                        continue
                    fi_v = float(fi_s); se_v = float(se_s)
                    if math.isnan(fi_v) or math.isnan(se_v):
                        continue
                    pt = [1.0 - fi_v, se_v]
                    if n_ < INIT_BUDGET:
                        Y_init_list.append(pt)
                    else:
                        Y_bo_list.append(pt)
                except (KeyError, ValueError):
                    continue
            Y_init_np = np.asarray(Y_init_list, dtype=float).reshape(-1, 2)
            Y_bo_np   = np.asarray(Y_bo_list,   dtype=float).reshape(-1, 2)
            latest, hist = plot_pareto(Y_init_np, Y_bo_np,
                                       hv_after, iter_idx, PLOT_DIR)
            print('  pareto plot -> %s  (history: %s)'
                  % (os.path.relpath(latest, BASE_DIR),
                     os.path.relpath(hist,   BASE_DIR)))
        except Exception as _e:
            print('  WARN: plot failed (%s)' % _e)

        print('  batch done.  HV: %.4f -> %.4f   (delta=%.4f)'
              % (hv_before, hv_after, hv_after - hv_before))
        bo_done += q_eff

    print('\n' + '=' * 64)
    print('  MBO loop exit.  BO specimens this run: %d' % bo_done)
    print('  See:  %s' % MBO_LOG_CSV)
    print('=' * 64)


if __name__ == '__main__':
    main()
