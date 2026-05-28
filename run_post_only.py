# -*- coding: utf-8 -*-
# =============================================================
#  run_post_only.py — Post-processing only (Abaqus 2024 compatible)
#
#  Use this when ODB files already exist (Phase 1 build + Phase 2
#  solver completed) and only Phase 3 (post-processing into
#  results.csv) needs to run.
#
#  Differences from run_all.py:
#    - Skips Phase 1 (build INP) and Phase 2 (solver submission).
#    - Wrapper script uses `exec(open(...).read())` instead of
#      Python-2-only `execfile()`, so it works on Abaqus 2024
#      (Python 3).
#    - Calls post_process_csv_ver3.py (Python 3 syntax).
#
#  Auto-discovers every `sim_<N>` folder under BASE_DIR that
#  contains `cutout_points_<N>.odb`. You can also restrict the
#  range via --start / --end.
#
#  Examples:
#    python run_post_only.py
#    python run_post_only.py --start 0 --end 11
#    python run_post_only.py --only 3,5,7
# =============================================================
import os
import re
import sys
import argparse
import subprocess
import time
import datetime

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
POST_SCRIPT = os.path.join(BASE_DIR, 'post_process_csv_ver3.py')


def write_temp_script(path, csv_num, target_script):
    """Wrapper that injects csv_num + BASE_DIR and exec()s the target.

    Abaqus 2024 ships Python 3, which removed `execfile`. Use
    exec(open(...).read(), globals()) instead so the target script
    runs in the wrapper's global scope (same effect as execfile).
    """
    with open(path, 'w') as f:
        f.write("# Auto-generated wrapper -- do not edit by hand.\n")
        f.write("csv_num = %d\n" % csv_num)
        f.write("BASE_DIR = r'%s'\n" % BASE_DIR)
        f.write("with open(r'%s', 'r') as _f:\n" % target_script)
        f.write("    exec(compile(_f.read(), r'%s', 'exec'), globals())\n"
                % target_script)


def run_abaqus_cae(script_path, label=""):
    cmd = ['abaqus', 'cae', 'noGUI=%s' % script_path]
    print("    CMD: %s" % ' '.join(cmd))
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        print("    WARNING: %s exited with code %d" % (label, rc))
    return rc


def discover_specimens(base_dir):
    """Return sorted list of N for every sim_<N>/cutout_points_<N>.odb."""
    pat = re.compile(r'^sim_(\d+)$')
    found = []
    for name in os.listdir(base_dir):
        m = pat.match(name)
        if not m:
            continue
        n = int(m.group(1))
        odb = os.path.join(base_dir, name, 'cutout_points_%d.odb' % n)
        if os.path.isfile(odb):
            found.append(n)
    found.sort()
    return found


def parse_args():
    ap = argparse.ArgumentParser(
        description='Post-process existing ODBs into results.csv. '
                    'Skips INP build and solver phases.')
    ap.add_argument('--start', type=int, default=None,
                    help='First specimen index (inclusive). Default: auto.')
    ap.add_argument('--end', type=int, default=None,
                    help='Last specimen index (inclusive). Default: auto.')
    ap.add_argument('--only', type=str, default=None,
                    help='Comma-separated specific indices, e.g. "3,5,7". '
                         'Overrides --start/--end.')
    ap.add_argument('--keep-temp', action='store_true',
                    help='Keep _temp_post_*.py wrappers for debugging.')
    return ap.parse_args()


def main():
    args = parse_args()

    discovered = discover_specimens(BASE_DIR)
    if not discovered:
        print("ERROR: No sim_<N>/cutout_points_<N>.odb found under %s"
              % BASE_DIR)
        sys.exit(2)

    if args.only:
        try:
            requested = sorted(int(s) for s in args.only.split(',') if s.strip())
        except ValueError:
            print("ERROR: --only must be comma-separated integers.")
            sys.exit(2)
        specimens = [n for n in requested if n in discovered]
        missing = [n for n in requested if n not in discovered]
        if missing:
            print("WARNING: requested but ODB missing: %s" % missing)
    else:
        lo = args.start if args.start is not None else discovered[0]
        hi = args.end   if args.end   is not None else discovered[-1]
        if hi < lo:
            print("ERROR: --end (%d) must be >= --start (%d)" % (hi, lo))
            sys.exit(2)
        specimens = [n for n in discovered if lo <= n <= hi]

    if not specimens:
        print("ERROR: No specimens to process after filtering.")
        sys.exit(2)

    if not os.path.isfile(POST_SCRIPT):
        print("ERROR: post-process script not found: %s" % POST_SCRIPT)
        sys.exit(2)

    total_start = time.time()
    print("=" * 60)
    print("  Post-Process Only (ODB -> results.csv)")
    print("  BASE_DIR  : %s" % BASE_DIR)
    print("  Specimens : %s  (%d total)" % (specimens, len(specimens)))
    print("  Script    : %s" % os.path.basename(POST_SCRIPT))
    print("  Start     : %s"
          % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 60)

    failures = []
    for n in specimens:
        print("\n  [Post] Specimen %d ..." % n)
        temp_file = os.path.join(BASE_DIR, '_temp_post_%d.py' % n)
        write_temp_script(temp_file, n, POST_SCRIPT)
        rc = run_abaqus_cae(temp_file, label="post_%d" % n)
        if rc != 0:
            failures.append(n)
        if not args.keep_temp:
            try:
                os.remove(temp_file)
            except OSError:
                pass
        print("  [Post] Specimen %d done (rc=%d)." % (n, rc))

    elapsed = time.time() - total_start
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)

    print("\n" + "=" * 60)
    print("  ALL DONE!")
    print("  Results : %s" % os.path.join(BASE_DIR, 'results.csv'))
    print("  Elapsed : %02d:%02d:%02d" % (h, m, s))
    if failures:
        print("  FAILED  : %s" % failures)
    else:
        print("  FAILED  : (none)")
    print("=" * 60)
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
