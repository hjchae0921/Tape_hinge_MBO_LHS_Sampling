# -*- coding: utf-8 -*-
"""
Thin wrappers around Abaqus CAE / Explicit invocations for the MBO loop.

Per-specimen pipeline (3 phases):
  Phase 1: tube_hinge_pinching.py     (build INP via abaqus cae noGUI)
  Phase 2: abaqus job=cutout_points_<N>  cpus=<C>  interactive
  Phase 3: post_process_csv_ver3.py   (ODB -> results.csv append/update)

q-batch execution (`run_batch_iteration`):
  Phase 1 (build)  : sequential   -- avoids Abaqus CAE license contention
                                     and races on shared inputs
  Phase 2 (solver) : parallel     -- q jobs concurrently, cpus_per_job each
                                     (default q=3, 4 cpus/job = 12 cores)
  Phase 3 (post)   : sequential   -- avoids racing writes on results.csv

`run_full_iteration` runs all 3 phases sequentially for ONE specimen and
is used as the retry path for any specimen that fails inside a batch.

Wrappers use exec(compile(..., 'exec'), globals()) instead of the
Python-2-only `execfile`, so they run on Abaqus 2024 (Python 3).
"""
import os
import subprocess

BUILD_SCRIPT = 'tube_hinge_pinching.py'
POST_SCRIPT  = 'post_process_csv_ver3.py'


def _wrapper_text(csv_num, base_dir, target_script):
    # `encoding='utf-8'` is required: target scripts may contain non-ASCII
    # characters (em-dash, Korean comments, ...), and Windows' default open()
    # decoder is cp949 in ko-KR locales -> UnicodeDecodeError on 0xE2 bytes.
    return (
        "# Auto-generated MBO wrapper -- do not edit by hand.\n"
        "csv_num = %d\n" % csv_num
        + "BASE_DIR = r'%s'\n" % base_dir
        + "with open(r'%s', 'r', encoding='utf-8') as _f:\n" % target_script
        + "    exec(compile(_f.read(), r'%s', 'exec'), globals())\n"
          % target_script
    )


def _write_wrapper(path, csv_num, base_dir, target_script):
    with open(path, 'w') as f:
        f.write(_wrapper_text(csv_num, base_dir, target_script))


def _run_cae(script_path):
    cmd = ['abaqus', 'cae', 'noGUI=%s' % script_path]
    print("    CMD:", ' '.join(cmd))
    return subprocess.call(cmd, shell=True)


def run_build_inp(csv_num, base_dir):
    target = os.path.join(base_dir, BUILD_SCRIPT)
    tmp = os.path.join(base_dir, '_mbo_phase1_%d.py' % csv_num)
    _write_wrapper(tmp, csv_num, base_dir, target)
    try:
        return _run_cae(tmp)
    finally:
        try: os.remove(tmp)
        except OSError: pass


def run_solver(csv_num, base_dir, cpus=4):
    sim_dir = os.path.join(base_dir, 'sim_%d' % csv_num)
    if not os.path.isdir(sim_dir):
        print('    ERROR: sim dir not found:', sim_dir)
        return 99
    job = 'cutout_points_%d' % csv_num
    cmd = ('abaqus job=%s cpus=%d double=both output_precision=full interactive'
           % (job, cpus))
    print("    CMD: %s   (cwd=%s)" % (cmd, sim_dir))
    return subprocess.call(cmd, cwd=sim_dir, shell=True)


def run_solver_parallel(csv_nums, base_dir, cpus_per_job=4):
    """
    Launch solver jobs for `csv_nums` in parallel and wait for all.

    Returns
    -------
    rcs : dict {csv_num: return_code}
          99 means the sim dir was missing (build did not run).
    """
    procs = []
    for n in csv_nums:
        sim_dir = os.path.join(base_dir, 'sim_%d' % n)
        if not os.path.isdir(sim_dir):
            print('    ERROR: sim dir not found for %d -> skipping launch' % n)
            procs.append((n, None, 99))
            continue
        job = 'cutout_points_%d' % n
        cmd = ('abaqus job=%s cpus=%d double=both output_precision=full interactive'
               % (job, cpus_per_job))
        print('    -> %s   (cwd=%s)' % (job, sim_dir))
        p = subprocess.Popen(cmd, cwd=sim_dir, shell=True)
        procs.append((n, p, None))

    rcs = {}
    for n, p, pre in procs:
        if p is None:
            rcs[n] = pre
            continue
        p.wait()
        rcs[n] = p.returncode
        print('    Specimen %d solver done (rc=%d)' % (n, p.returncode))
    return rcs


def run_post(csv_num, base_dir):
    target = os.path.join(base_dir, POST_SCRIPT)
    tmp = os.path.join(base_dir, '_mbo_phase3_%d.py' % csv_num)
    _write_wrapper(tmp, csv_num, base_dir, target)
    try:
        return _run_cae(tmp)
    finally:
        try: os.remove(tmp)
        except OSError: pass


def run_full_iteration(csv_num, base_dir, cpus=4):
    """
    Run phases 1->2->3 for ONE specimen sequentially. Used for retries.

    Returns
    -------
    success : bool
    phase   : str   ('build' | 'solver' | 'post' | 'done')
    rcs     : dict  return-codes per phase
    """
    rcs = {}

    rcs['build'] = run_build_inp(csv_num, base_dir)
    if rcs['build'] != 0:
        return False, 'build', rcs

    rcs['solver'] = run_solver(csv_num, base_dir, cpus=cpus)
    if rcs['solver'] != 0:
        return False, 'solver', rcs

    odb = os.path.join(base_dir, 'sim_%d' % csv_num,
                       'cutout_points_%d.odb' % csv_num)
    if not os.path.isfile(odb):
        rcs['odb_missing'] = True
        return False, 'solver', rcs

    rcs['post'] = run_post(csv_num, base_dir)
    if rcs['post'] != 0:
        return False, 'post', rcs

    return True, 'done', rcs


def run_batch_iteration(csv_nums, base_dir, cpus_per_job=4):
    """
    Run a q-batch through the 3-phase pipeline.

    Phase 1 (build INP)        : sequential (abaqus cae license + repo writes)
    Phase 2 (Abaqus Explicit)  : parallel  (cpus_per_job each, len(csv_nums) jobs)
    Phase 3 (post-process ODB) : sequential (writes shared results.csv)

    Returns
    -------
    statuses : dict {csv_num: (success: bool, last_phase: str, rcs: dict)}
    """
    statuses = {}

    # ---- Phase 1: build INP (sequential) ----
    print('  [batch phase1] build INP for %s' % csv_nums)
    for n in csv_nums:
        rc = run_build_inp(n, base_dir)
        if rc != 0:
            statuses[n] = (False, 'build', {'build': rc})
        else:
            statuses[n] = (None, None, {'build': 0})

    # ---- Phase 2: solver (parallel) ----
    pending = [n for n in csv_nums if statuses[n][0] is None]
    print('  [batch phase2] solver parallel for %s   cpus/job=%d'
          % (pending, cpus_per_job))
    if pending:
        rcs_solver = run_solver_parallel(pending, base_dir,
                                         cpus_per_job=cpus_per_job)
        for n in pending:
            rc_s = rcs_solver[n]
            statuses[n][2]['solver'] = rc_s
            if rc_s != 0:
                statuses[n] = (False, 'solver', statuses[n][2])
                continue
            odb = os.path.join(base_dir, 'sim_%d' % n,
                               'cutout_points_%d.odb' % n)
            if not os.path.isfile(odb):
                statuses[n][2]['odb_missing'] = True
                statuses[n] = (False, 'solver', statuses[n][2])

    # ---- Phase 3: post-process (sequential) ----
    pending = [n for n in csv_nums if statuses[n][0] is None]
    print('  [batch phase3] post-process for %s' % pending)
    for n in pending:
        rc_p = run_post(n, base_dir)
        statuses[n][2]['post'] = rc_p
        if rc_p != 0:
            statuses[n] = (False, 'post', statuses[n][2])
        else:
            statuses[n] = (True, 'done', statuses[n][2])

    return statuses
