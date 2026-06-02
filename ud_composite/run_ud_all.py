# -*- coding: utf-8 -*-
# =============================================================
#  run_ud_all.py — UD composite batch driver
#
#  For each layup in --layup-start..--layup-end :
#    Phase 1: Build INP for specimen 291, then 197         (sequential CAE)
#    Phase 2: Run 291 + 197 solvers in parallel (6 cpus each, Explicit)
#    Phase 3: Post-process ODBs sequentially → results_ud.csv
#
#  Targets Abaqus 2017 (Python 2.7) — uses execfile() in temp wrappers.
#
#  Examples (PC = 12 cores):
#    abaqus python run_ud_all.py                          # all 74 layups
#    abaqus python run_ud_all.py --layup-start 36 --layup-end 54
#    abaqus python run_ud_all.py --layup-ids 11,36,55,69
#    abaqus python run_ud_all.py --specimens 291         # only one specimen
#
#  Layup IDs (see ud_layups.py / ud_ply_angle_candidates.md):
#    1..20  : 2-ply,  21..35 : 3-ply,  36..54 : 4-ply,
#    55..68 : 6-ply,  69..74 : 8-ply.
# =============================================================
from __future__ import print_function
import os
import sys
import argparse
import subprocess
import time
import datetime

# =============================================================
# Paths
# =============================================================
UD_DIR        = os.path.dirname(os.path.abspath(__file__))
BASE_DIR      = os.path.dirname(UD_DIR)
BUILD_SCRIPT  = os.path.join(UD_DIR, 'build_ud_specimen.py')
POST_SCRIPT   = os.path.join(UD_DIR, 'post_ud_specimen.py')

if UD_DIR not in sys.path:
    sys.path.insert(0, UD_DIR)
from ud_layups import LAYUPS, get_layup


# =============================================================
# Helpers
# =============================================================
def write_temp_script(path, csv_num, layup_id, target_script):
    """Wrapper that injects csv_num + layup_id + BASE_DIR then execfiles target."""
    f = open(path, 'w')
    try:
        f.write("csv_num   = %d\n" % csv_num)
        f.write("layup_id  = %d\n" % layup_id)
        f.write("BASE_DIR  = r'%s'\n" % BASE_DIR)
        f.write("execfile(r'%s')\n" % target_script)
    finally:
        f.close()


def run_abaqus_cae(script_path, label=''):
    cmd = ['abaqus', 'cae', 'noGUI=%s' % script_path]
    print('    CMD: %s' % ' '.join(cmd))
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        print('    WARNING: %s exited rc=%d' % (label, rc))
    return rc


def run_solvers_parallel(jobs, cpus_per_job):
    """jobs = [(csv_num, layup_id, sim_dir, job_name), ...]"""
    procs = []
    for (csv_num, layup_id, sim_dir, job_name) in jobs:
        cmd = ('abaqus job=%s cpus=%d double=both '
               'output_precision=full interactive'
               % (job_name, cpus_per_job))
        print('    -> %s   (dir: %s)' % (job_name, sim_dir))
        p = subprocess.Popen(cmd, cwd=sim_dir, shell=True)
        procs.append((csv_num, layup_id, p))

    for (csv_num, layup_id, p) in procs:
        p.wait()
        print('    L%d / S%d solver done (rc=%d)'
              % (layup_id, csv_num, p.returncode))


def parse_ids(s):
    out = []
    for tok in s.split(','):
        tok = tok.strip()
        if tok:
            out.append(int(tok))
    return out


# =============================================================
# Main
# =============================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description='UD composite batch automation (74 layups x 2 specimens).')
    ap.add_argument('--layup-start', type=int, default=1,
                    help='First layup ID (default: 1)')
    ap.add_argument('--layup-end', type=int, default=len(LAYUPS),
                    help='Last layup ID, inclusive (default: %d)' % len(LAYUPS))
    ap.add_argument('--layup-ids', type=str, default='',
                    help='Explicit comma-separated layup IDs (overrides start/end).')
    ap.add_argument('--specimens', type=str, default='291,197',
                    help='Comma-separated specimen numbers (default: 291,197).')
    ap.add_argument('--cpus', type=int, default=6,
                    help='CPUs per Abaqus job (default: 6).')
    return ap.parse_args()


def main():
    args = parse_args()

    if args.layup_ids.strip():
        layup_ids = parse_ids(args.layup_ids)
    else:
        if args.layup_end < args.layup_start:
            print('ERROR: --layup-end (%d) must be >= --layup-start (%d)'
                  % (args.layup_end, args.layup_start))
            sys.exit(2)
        layup_ids = list(range(args.layup_start, args.layup_end + 1))

    specimens = parse_ids(args.specimens)
    if not specimens:
        print('ERROR: --specimens parsed empty.')
        sys.exit(2)

    total_start = time.time()
    print('=' * 60)
    print('  UD Composite Batch Automation')
    print('  Layups   : %d entries  -> %s' % (len(layup_ids), layup_ids))
    print('  Specimens: %s' % specimens)
    print('  CPUs/job : %d   (2 parallel per layup -> %d cores total)'
          % (args.cpus, args.cpus * len(specimens)))
    print('  Start    : %s'
          % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    print('=' * 60)

    for batch_idx, layup_id in enumerate(layup_ids):
        lay = get_layup(layup_id)
        print('')
        print('=' * 60)
        print('  LAYUP %d / %d   id=%d  %s   %s   angles=%s'
              % (batch_idx + 1, len(layup_ids),
                 layup_id, lay['group'], lay['label'], lay['angles']))
        print('=' * 60)

        # -------------------------------------------------------
        # Phase 1: Build INP files for both specimens (sequential)
        # -------------------------------------------------------
        for csv_num in specimens:
            print('  [Phase1] Build INP  L%d / S%d ...' % (layup_id, csv_num))
            temp_file = os.path.join(
                UD_DIR, '_temp_build_L%d_S%d.py' % (layup_id, csv_num))
            write_temp_script(temp_file, csv_num, layup_id, BUILD_SCRIPT)
            run_abaqus_cae(temp_file, label='build_L%d_S%d' % (layup_id, csv_num))
            try:
                os.remove(temp_file)
            except OSError:
                pass
            print('  [Phase1] INP done   L%d / S%d' % (layup_id, csv_num))

        # -------------------------------------------------------
        # Phase 2: Submit solver jobs in parallel
        # -------------------------------------------------------
        jobs = []
        for csv_num in specimens:
            sim_dir  = os.path.join(UD_DIR,
                                    'sim_L%d_S%d' % (layup_id, csv_num))
            job_name = 'ud_L%d_S%d' % (layup_id, csv_num)
            jobs.append((csv_num, layup_id, sim_dir, job_name))

        print('  [Phase2] Launching %d solver(s) in parallel (cpus=%d each) ...'
              % (len(jobs), args.cpus))
        run_solvers_parallel(jobs, args.cpus)
        print('  [Phase2] All solvers for layup %d complete.' % layup_id)

        # -------------------------------------------------------
        # Phase 3: Post-process ODBs (sequential, results_ud.csv writer
        # is non-reentrant)
        # -------------------------------------------------------
        for csv_num in specimens:
            print('  [Phase3] Post-process L%d / S%d ...'
                  % (layup_id, csv_num))
            temp_file = os.path.join(
                UD_DIR, '_temp_post_L%d_S%d.py' % (layup_id, csv_num))
            write_temp_script(temp_file, csv_num, layup_id, POST_SCRIPT)
            run_abaqus_cae(temp_file, label='post_L%d_S%d' % (layup_id, csv_num))
            try:
                os.remove(temp_file)
            except OSError:
                pass
            print('  [Phase3] Done       L%d / S%d' % (layup_id, csv_num))

        print('  LAYUP %d COMPLETE.' % layup_id)

    elapsed = time.time() - total_start
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    print('')
    print('=' * 60)
    print('  ALL DONE!')
    print('  Results : %s'
          % os.path.join(UD_DIR, 'results_ud.csv'))
    print('  Elapsed : %02d:%02d:%02d' % (h, m, s))
    print('=' * 60)


if __name__ == '__main__':
    main()
