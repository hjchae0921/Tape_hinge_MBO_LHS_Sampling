# -*- coding: utf-8 -*-
# =============================================================
#  run_all.py — Pinching + Folding batch automation
#
#  Per-machine CLI arguments (configure for each PC):
#    --start  N  : first specimen index (default 0)
#    --end    N  : last  specimen index (default 39, inclusive)
#    --batch  N  : concurrent solver jobs per batch (default 2)
#    --cpus   N  : CPUs per Abaqus job (default 6)
#
#  Flow per batch:
#    Phase 1: Build INP files (sequential, Abaqus CAE)
#    Phase 2: Submit solver jobs in parallel (Abaqus Explicit)
#    Phase 3: Post-process ODBs (sequential, Abaqus CAE)
#
#  Examples:
#    python run_all.py --start 0  --end 18 --batch 3        # PC1 (6c x 3job)
#    python run_all.py --start 19 --end 30 --batch 2        # PC2 (6c x 2job)
#    python run_all.py --start 31 --end 39 --batch 2        # PC3 (6c x 2job)
# =============================================================
import os
import sys
import argparse
import subprocess
import time
import datetime

# =============================================================
# Paths (script-relative)
# =============================================================
BASE_DIR     = os.path.dirname(os.path.abspath(__file__))
BUILD_SCRIPT = os.path.join(BASE_DIR, 'tube_hinge_pinching.py')
POST_SCRIPT  = os.path.join(BASE_DIR, 'post_process_csv_ver2.py')


def write_temp_script(path, csv_num, target_script):
    """Generate a wrapper that injects csv_num + BASE_DIR before execfiling target."""
    with open(path, 'w') as f:
        f.write("csv_num = %d\n" % csv_num)
        f.write("BASE_DIR = r'%s'\n" % BASE_DIR)
        f.write("execfile(r'%s')\n" % target_script)


def run_abaqus_cae(script_path, label=""):
    """Run abaqus cae noGUI=<script> and wait for completion."""
    cmd = ['abaqus', 'cae', 'noGUI=%s' % script_path]
    print("    CMD: %s" % ' '.join(cmd))
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        print("    WARNING: %s exited with code %d" % (label, rc))
    return rc


def run_solvers_parallel(batch_nums, cpus_per_job):
    """Launch solver jobs in parallel and wait for all to finish."""
    procs = []
    for n in batch_nums:
        sim_dir  = os.path.join(BASE_DIR, 'sim_%d' % n)
        job_name = 'cutout_points_%d' % n
        cmd = ('abaqus job=%s cpus=%d double=both output_precision=full interactive'
               % (job_name, cpus_per_job))
        print("    -> %s  (dir: %s)" % (job_name, sim_dir))
        p = subprocess.Popen(cmd, cwd=sim_dir, shell=True)
        procs.append((n, p))

    for n, p in procs:
        p.wait()
        print("    Specimen %d solver done (rc=%d)" % (n, p.returncode))


# =============================================================
# Main
# =============================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description='Pinching+Folding batch automation. Set per-machine '
                    'specimen range and parallelism via CLI.')
    ap.add_argument('--start', type=int, default=0,
                    help='First specimen index (default: 0)')
    ap.add_argument('--end', type=int, default=39,
                    help='Last specimen index, inclusive (default: 39)')
    ap.add_argument('--batch', type=int, default=2,
                    help='Concurrent solver jobs per batch (default: 2)')
    ap.add_argument('--cpus', type=int, default=6,
                    help='CPUs per Abaqus job (default: 6)')
    return ap.parse_args()


def main():
    args = parse_args()
    if args.end < args.start:
        print("ERROR: --end (%d) must be >= --start (%d)" % (args.end, args.start))
        sys.exit(2)

    total_start = time.time()
    specimens = list(range(args.start, args.end + 1))
    total_batches = (len(specimens) + args.batch - 1) // args.batch

    print("=" * 60)
    print("  Pinching Batch Automation")
    print("  Specimens : %d .. %d (%d total)" % (args.start, args.end, len(specimens)))
    print("  Batch size: %d   CPUs/job: %d" % (args.batch, args.cpus))
    print("  Start     : %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 60)

    for batch_idx in range(total_batches):
        i_start = batch_idx * args.batch
        batch_nums = specimens[i_start : i_start + args.batch]

        print("\n" + "=" * 60)
        print("  BATCH %d / %d  (specimens %s)"
              % (batch_idx + 1, total_batches, batch_nums))
        print("=" * 60)

        # -----------------------------------------------------
        # Phase 1: Build INP files (sequential)
        # -----------------------------------------------------
        for n in batch_nums:
            print("  [Phase1] Building INP for specimen %d ..." % n)
            temp_file = os.path.join(BASE_DIR, '_temp_build_%d.py' % n)
            write_temp_script(temp_file, n, BUILD_SCRIPT)
            run_abaqus_cae(temp_file, label="build_%d" % n)
            try:
                os.remove(temp_file)
            except OSError:
                pass
            print("  [Phase1] Specimen %d INP done." % n)

        # -----------------------------------------------------
        # Phase 2: Submit solver jobs in parallel
        # -----------------------------------------------------
        print("  [Phase2] Submitting %d solver jobs in parallel ..."
              % len(batch_nums))
        run_solvers_parallel(batch_nums, args.cpus)
        print("  [Phase2] All solvers in batch %d complete." % (batch_idx + 1))

        # -----------------------------------------------------
        # Phase 3: Post-process ODBs (sequential)
        # -----------------------------------------------------
        for n in batch_nums:
            print("  [Phase3] Post-processing specimen %d ..." % n)
            temp_file = os.path.join(BASE_DIR, '_temp_post_%d.py' % n)
            write_temp_script(temp_file, n, POST_SCRIPT)
            run_abaqus_cae(temp_file, label="post_%d" % n)
            try:
                os.remove(temp_file)
            except OSError:
                pass
            print("  [Phase3] Specimen %d post-processed." % n)

        print("  BATCH %d COMPLETE" % (batch_idx + 1))

    elapsed = time.time() - total_start
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)

    print("\n" + "=" * 60)
    print("  ALL DONE!")
    print("  Results : %s" % os.path.join(BASE_DIR, 'results.csv'))
    print("  Elapsed : %02d:%02d:%02d" % (h, m, s))
    print("=" * 60)


if __name__ == '__main__':
    main()
