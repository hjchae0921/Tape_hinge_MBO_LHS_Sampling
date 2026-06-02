# -*- coding: utf-8 -*-
# =============================================================
#  post_ud_specimen.py — UD CompositeShell post-processor
#
#  Opens ud_composite/sim_L<layup_id>_S<csv_num>/ud_L*_S*.odb,
#  scans every frame of (pinching, folding) over every element /
#  integration point / section point, and collects the global max
#  of the built-in CFAILURE components:
#
#      TSAIH  -- Tsai-Hill
#      TSAIW  -- Tsai-Wu
#      MSTRS  -- Maximum Stress
#
#  Results are appended/updated to ud_composite/results_ud.csv,
#  keyed on LAYUP_ID, with one row per layup_id containing both
#  specimens (291 & 197).
#
#  Runner injects:
#      csv_num   : 291 or 197
#      layup_id  : 1..74
#      BASE_DIR  : repo root
#
#  Run via runner — `abaqus cae noGUI=_temp_post_L<id>_S<num>.py`.
#  Compatible with Abaqus 2017 (Python 2.7) and 2024 (Python 3).
# =============================================================
from __future__ import print_function
import os
import sys
import csv

from abaqus import session

# =============================================================
# 0. Inputs from runner
# =============================================================
try:
    csv_num
except NameError:
    csv_num = 291
try:
    layup_id
except NameError:
    layup_id = 11

try:
    BASE_DIR
except NameError:
    try:
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        BASE_DIR = os.getcwd()

UD_DIR      = os.path.join(BASE_DIR, 'ud_composite')
SIM_DIR     = os.path.join(UD_DIR, 'sim_L%d_S%d' % (layup_id, csv_num))
RESULTS_CSV = os.path.join(UD_DIR, 'results_ud.csv')

odb_path = os.path.join(SIM_DIR, 'ud_L%d_S%d.odb' % (layup_id, csv_num))

if UD_DIR not in sys.path:
    sys.path.insert(0, UD_DIR)
from ud_layups import get_layup, PLY_THICKNESS

LAY = get_layup(layup_id)
PLY_ANGLES = LAY['angles']
N_PLIES    = len(PLY_ANGLES)
TOT_THICK  = N_PLIES * PLY_THICKNESS

print('============================================')
print('  UD POST  layup_id=%d (%s)  specimen=%d' % (layup_id, LAY['label'], csv_num))
print('  ODB: %s' % odb_path)
print('============================================')

# =============================================================
# 1. Open ODB
# =============================================================
if odb_path in session.odbs.keys():
    session.odbs[odb_path].close()
odb = session.openOdb(name=odb_path, readOnly=True)

inst_keys = [k for k in odb.rootAssembly.instances.keys() if k.upper() != 'ASSEMBLY']
if not inst_keys:
    raise KeyError('No part instance found in ODB.')
inst_obj = odb.rootAssembly.instances[inst_keys[0]]
print('Using instance: %s' % inst_keys[0])

# =============================================================
# 2. Scan CFAILURE across all frames / steps
# =============================================================
TARGET_STEPS = ['pinching', 'folding']
available_steps = list(odb.steps.keys())
step_names = [s for s in TARGET_STEPS if s in available_steps]
print('Processing steps: %s' % str(step_names))


def _find_indices(labels):
    """Return (idx_TH, idx_TW, idx_MS) into a CFAILURE component-label tuple."""
    name_map = {}
    for i, lbl in enumerate(labels):
        name_map[lbl.upper()] = i
    return (
        name_map.get('TSAIH'),
        name_map.get('TSAIW'),
        name_map.get('MSTRS'),
    )


def _safe_max(current, candidate):
    if candidate is None:
        return current
    if current is None or candidate > current:
        return candidate
    return current


max_TH = None
max_TW = None
max_MS = None
n_values = 0

for step_name in step_names:
    step = odb.steps[step_name]
    nFrames = len(step.frames)
    print('--- Step %s : %d frames ---' % (step_name, nFrames))

    step_TH = None
    step_TW = None
    step_MS = None

    for frame_index, frame in enumerate(step.frames):
        if 'CFAILURE' not in frame.fieldOutputs.keys():
            continue
        field = frame.fieldOutputs['CFAILURE']
        sub   = field.getSubset(region=inst_obj)
        labels = sub.componentLabels
        idx_TH, idx_TW, idx_MS = _find_indices(labels)

        for v in sub.values:
            d = v.data
            if idx_TH is not None:
                step_TH = _safe_max(step_TH, d[idx_TH])
            if idx_TW is not None:
                step_TW = _safe_max(step_TW, d[idx_TW])
            if idx_MS is not None:
                step_MS = _safe_max(step_MS, d[idx_MS])
            n_values += 1

        if (frame_index + 1) % 50 == 0:
            print('  frame %d/%d  TH=%s TW=%s MS=%s'
                  % (frame_index + 1, nFrames,
                     step_TH, step_TW, step_MS))

    print('  Step %s: TH=%s  TW=%s  MS=%s'
          % (step_name, step_TH, step_TW, step_MS))

    max_TH = _safe_max(max_TH, step_TH)
    max_TW = _safe_max(max_TW, step_TW)
    max_MS = _safe_max(max_MS, step_MS)

print('')
print('Scanned %d CFAILURE values.' % n_values)
print('GLOBAL  Tsai-Hill = %s' % max_TH)
print('GLOBAL  Tsai-Wu   = %s' % max_TW)
print('GLOBAL  Max stress= %s' % max_MS)

session.odbs[odb_path].close()

# =============================================================
# 3. Update results_ud.csv (one row per layup_id)
# =============================================================
def _fmt(v):
    return '' if v is None else float(v)


HEADER = [
    'LAYUP_ID', 'GROUP', 'LAYUP_LABEL', 'PLY_COUNT', 'LAYUP_ANGLES',
    'TSAI_HILL_291', 'TSAI_HILL_197',
    'TSAI_WU_291',   'TSAI_WU_197',
    'MAX_STRESS_291','MAX_STRESS_197',
]

th_col = 'TSAI_HILL_%d'   % csv_num
tw_col = 'TSAI_WU_%d'     % csv_num
ms_col = 'MAX_STRESS_%d'  % csv_num

# Read existing rows
existing = []
if os.path.exists(RESULTS_CSV):
    try:
        f = open(RESULTS_CSV, 'r')
        try:
            reader = csv.DictReader(f)
            for row in reader:
                existing.append(row)
        finally:
            f.close()
    except Exception as e:
        print('WARNING: could not read existing %s: %s' % (RESULTS_CSV, str(e)))

# Find or create the row for this layup
target_row = None
for row in existing:
    if str(row.get('LAYUP_ID', '')).strip() == str(layup_id):
        target_row = row
        break

if target_row is None:
    target_row = {}
    for h in HEADER:
        target_row[h] = ''
    target_row['LAYUP_ID']      = layup_id
    target_row['GROUP']         = LAY['group']
    target_row['LAYUP_LABEL']   = LAY['label']
    target_row['PLY_COUNT']     = N_PLIES
    target_row['LAYUP_ANGLES']  = '/'.join(str(a) for a in PLY_ANGLES)
    existing.append(target_row)
else:
    target_row['GROUP']        = LAY['group']
    target_row['LAYUP_LABEL']  = LAY['label']
    target_row['PLY_COUNT']    = N_PLIES
    target_row['LAYUP_ANGLES'] = '/'.join(str(a) for a in PLY_ANGLES)

target_row[th_col] = _fmt(max_TH)
target_row[tw_col] = _fmt(max_TW)
target_row[ms_col] = _fmt(max_MS)

# Sort by LAYUP_ID (numeric, fall back to string)
def _sort_key(r):
    try:
        return (0, int(r.get('LAYUP_ID', 0)))
    except Exception:
        return (1, str(r.get('LAYUP_ID', '')))
existing.sort(key=_sort_key)

# Python 2 vs 3 csv open mode
try:
    f = open(RESULTS_CSV, 'w', newline='')
except TypeError:
    f = open(RESULTS_CSV, 'wb')
try:
    writer = csv.DictWriter(f, fieldnames=HEADER, lineterminator='\n',
                            extrasaction='ignore')
    writer.writeheader()
    writer.writerows(existing)
finally:
    f.close()

print('')
print('============================================')
print('  %s updated' % RESULTS_CSV)
print('  LAYUP_ID=%d  specimen=%d' % (layup_id, csv_num))
print('    %s = %s' % (th_col, target_row[th_col]))
print('    %s = %s' % (tw_col, target_row[tw_col]))
print('    %s = %s' % (ms_col, target_row[ms_col]))
print('============================================')
