# -*- coding: utf-8 -*-
# =============================================================
#  post_process_csv_ver3.py  (Abaqus 2024 / Python 3 compatible)
#
#  ver2 -> ver3 changes:
#    - All `print '...'` statements converted to print() function calls
#      (Abaqus 2024 ships Python 3; the old syntax is a SyntaxError).
#    - DictWriter file opened with newline='' to avoid extra blank lines
#      that Python 3 csv module produces on Windows.
#    - No logic change in the FI / SE computation.
#
#  Run:
#    abaqus cae noGUI=post_process_csv_ver3.py -- csv_num=<N>
#  or via wrapper that sets csv_num + BASE_DIR before exec().
# =============================================================
from __future__ import print_function
import os
import csv
import math

from abaqus import mdb, session
from abaqusConstants import *
from math import sqrt as _sqrt

# =============================================================
# User input
# =============================================================
try:
    csv_num
except NameError:
    csv_num = 0

try:
    BASE_DIR
except NameError:
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        BASE_DIR = os.getcwd()

SPECIMEN_DIR = os.path.join(BASE_DIR, 'specimen')
SIM_DIR      = os.path.join(BASE_DIR, 'sim_%d' % csv_num)
RESULTS_CSV  = os.path.join(BASE_DIR, 'results.csv')

odbPath      = os.path.join(SIM_DIR, 'cutout_points_%d.odb' % csv_num)
csv_specimen = os.path.join(SPECIMEN_DIR, 'cutout_points_%d.csv' % csv_num)

print('============================================')
print('  Post-processing specimen %d' % csv_num)
print('  ODB: %s' % odbPath)
print('============================================')

# =============================================================
# 1. Read design variables from cutout CSV
# =============================================================
eccentricity = ''
y1 = ''; x1 = ''; Dx = ''; y3 = ''; yP2 = ''; fP3 = ''
x2_derived = ''; x3_derived = ''; xP3_derived = ''; yP3_derived = ''

try:
    with open(csv_specimen, 'r') as f:
        reader = csv.reader(f)
        rows = list(reader)
    if len(rows) > 1:
        data_row = rows[1]
        if len(data_row) >= 3 and data_row[2].strip():
            eccentricity = float(data_row[2].strip())
        if len(data_row) >= 9:
            if data_row[3].strip(): y1  = float(data_row[3].strip())
            if data_row[4].strip(): x1  = float(data_row[4].strip())
            if data_row[5].strip(): Dx  = float(data_row[5].strip())
            if data_row[6].strip(): y3  = float(data_row[6].strip())
            if data_row[7].strip(): yP2 = float(data_row[7].strip())
            if data_row[8].strip(): fP3 = float(data_row[8].strip())
        if len(data_row) >= 13:
            if data_row[9].strip():  x2_derived  = float(data_row[9].strip())
            if data_row[10].strip(): x3_derived  = float(data_row[10].strip())
            if data_row[11].strip(): xP3_derived = float(data_row[11].strip())
            if data_row[12].strip(): yP3_derived = float(data_row[12].strip())
    print('ECCENTRICITY: %s' % eccentricity)
    print('DVs: y1=%s x1=%s Dx=%s y3=%s yP2=%s fP3=%s' % (y1, x1, Dx, y3, yP2, fP3))
    print('Derived: x2=%s x3=%s xP3=%s yP3=%s' % (x2_derived, x3_derived, xP3_derived, yP3_derived))
except Exception as e:
    print('WARNING: Failed to read cutout CSV - %s' % str(e))

# =============================================================
# 2. Open ODB
# =============================================================
if odbPath in session.odbs.keys():
    session.odbs[odbPath].close()

odb = session.openOdb(name=odbPath, readOnly=True)
print('ODB opened: %s' % odbPath)

inst_keys = [k for k in odb.rootAssembly.instances.keys()
             if k.upper() != 'ASSEMBLY']
if len(inst_keys) == 0:
    raise KeyError('No part instance found in ODB.')
inst_obj = odb.rootAssembly.instances[inst_keys[0]]
print('Using instance: %s' % inst_keys[0])

# =============================================================
# 3. Failure strength parameters
# =============================================================
sF1t = 139.47;  sF1c = 63.42
sF3  = 17.73;   sF4  = 5.07;   sF6 = 1.53
sJ1  = 1.0/sF1t - 1.0/sF1c
sK11 = 1.0/(sF1t*sF1c)
sK33 = 1.0/(sF3**2)
sK44 = 1.0/(sF4**2)
sK66 = 1.0/(sF6**2)
sK12 = -sK11/2.0


def _select_root(N_self, root_plus, root_minus):
    if N_self >= 0:
        candidates = [r for r in (root_plus, root_minus) if r > 0]
    else:
        candidates = [r for r in (root_plus, root_minus) if r < 0]

    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return candidates[0] if abs(candidates[0]) >= abs(candidates[1]) else candidates[1]


# =============================================================
# 4. Compute MAX_FI per step (pinching / folding / all)
# =============================================================
available_steps = list(odb.steps.keys())
target_steps = ['pinching', 'folding']
stepNames = [s for s in target_steps if s in available_steps]
print('Processing steps: %s' % str(stepNames))

max_fi_all      = -1e99
max_fi_pinching = -1e99
max_fi_folding  = -1e99

for step_name in stepNames:
    step = odb.steps[step_name]
    nFrames = len(step.frames)
    print('--- Step: %s (%d frames) ---' % (step_name, nFrames))

    step_max = -1e99

    for frame_index, frame in enumerate(step.frames):
        available = list(frame.fieldOutputs.keys())
        if 'SF' not in available or 'SM' not in available:
            continue

        fValues = frame.fieldOutputs['SF'].getSubset(region=inst_obj).values
        mValues = frame.fieldOutputs['SM'].getSubset(region=inst_obj).values

        mMap = {}
        for mv in mValues:
            key = (int(mv.elementLabel), int(mv.integrationPoint))
            mMap[key] = mv.data

        for fv in fValues:
            key = (int(fv.elementLabel), int(fv.integrationPoint))
            if key not in mMap:
                continue
            N1, N2, N3, N12, N13, N23 = fv.data
            M1, M2, M12               = mMap[key]

            Nx  =  0.5*(N1 + N2) + N12
            Ny  =  0.5*(N1 + N2) - N12
            Nxy =  0.5*(N2 - N1)
            Mx  =  0.5*(M1 + M2) + M12
            My  =  0.5*(M1 + M2) - M12
            Mxy =  0.5*(M2 - M1)

            M = Mx if abs(Mx) >= abs(My) else My

            fIndexIP = sJ1*(Nx+Ny) + sK11*(Nx**2 + Ny**2) + sK12*Nx*Ny + sK33*Nxy**2

            if fIndexIP >= 1.0:
                fIndexC = 1.0 + abs(M)/sF4
            else:
                dx = (sJ1 + sK12*Ny)**2 - 4*sK11*(sJ1*Ny + sK11*Ny**2 + sK33*Nxy**2 - 1)
                if dx >= 0:
                    sqrt_dx = _sqrt(dx)
                    Fx_plus  = (-(sJ1+sK12*Ny) + sqrt_dx)/(2*sK11)
                    Fx_minus = (-(sJ1+sK12*Ny) - sqrt_dx)/(2*sK11)
                    sFx = _select_root(Nx, Fx_plus, Fx_minus)
                else:
                    sFx = None

                if sFx is None or abs(sFx) < 1e-9:
                    fIndexCx = fIndexIP + abs(M)/sF4
                else:
                    fIndexCx = abs(Nx/sFx) + abs(M)/sF4

                dy = (sJ1 + sK12*Nx)**2 - 4*sK11*(sJ1*Nx + sK11*Nx**2 + sK33*Nxy**2 - 1)
                if dy >= 0:
                    sqrt_dy = _sqrt(dy)
                    Fy_plus  = (-(sJ1+sK12*Nx) + sqrt_dy)/(2*sK11)
                    Fy_minus = (-(sJ1+sK12*Nx) - sqrt_dy)/(2*sK11)
                    sFy = _select_root(Ny, Fy_plus, Fy_minus)
                else:
                    sFy = None

                if sFy is None or abs(sFy) < 1e-9:
                    fIndexCy = fIndexIP + abs(M)/sF4
                else:
                    fIndexCy = abs(Ny/sFy) + abs(M)/sF4

                fIndexC = fIndexCx if fIndexCx >= fIndexCy else fIndexCy

            if fIndexC > step_max:
                step_max = fIndexC

        if (frame_index + 1) % 50 == 0:
            print('  frame %d/%d  step_max=%.4f' % (frame_index + 1, nFrames, step_max))

    if step_name == 'pinching' and step_max > -1e99:
        max_fi_pinching = step_max
    if step_name == 'folding' and step_max > -1e99:
        max_fi_folding = step_max
    if step_max > max_fi_all:
        max_fi_all = step_max

    print('  Step %s MAX_FI = %.6f' % (step_name, step_max if step_max > -1e99 else 0))

if max_fi_all      <= -1e99: max_fi_all      = ''
if max_fi_pinching <= -1e99: max_fi_pinching = ''
if max_fi_folding  <= -1e99: max_fi_folding  = ''

print('')
print('MAX_FI_ALL      : %s' % max_fi_all)
print('MAX_FI_PINCHING : %s' % max_fi_pinching)
print('MAX_FI_FOLDING  : %s' % max_fi_folding)

# =============================================================
# 5. MAX_SE -- max ALLSE across all steps (History Output)
# =============================================================
max_se = ''
try:
    all_se_values = []
    for step_name in stepNames:
        step = odb.steps[step_name]
        for region_name in step.historyRegions.keys():
            region = step.historyRegions[region_name]
            if 'ALLSE' in region.historyOutputs:
                ho = region.historyOutputs['ALLSE']
                for time_val, data_val in ho.data:
                    all_se_values.append(data_val)
    if all_se_values:
        max_se = max(all_se_values)
        print('MAX_SE (ALLSE) : %s' % max_se)
    else:
        print('WARNING: ALLSE history data not found.')
except Exception as e:
    print('WARNING: Failed to extract MAX_SE - %s' % str(e))

# =============================================================
# 6. Close ODB
# =============================================================
session.odbs[odbPath].close()
print('ODB closed.')

# =============================================================
# 7. Write / update results.csv
# =============================================================
HEADER = [
    'SPECIMEN_NUM',
    'MAX_FI_ALL', 'MAX_FI_PINCHING', 'MAX_FI_FOLDING',
    'MAX_SE',
    'Eccentricity', 'y1', 'x1', 'Dx', 'y3', 'yP2', 'fP3',
    'x2_derived', 'x3_derived', 'xP3_derived', 'yP3_derived',
]

new_row = {
    'SPECIMEN_NUM'    : csv_num,
    'MAX_FI_ALL'      : max_fi_all,
    'MAX_FI_PINCHING' : max_fi_pinching,
    'MAX_FI_FOLDING'  : max_fi_folding,
    'MAX_SE'          : max_se,
    'Eccentricity'    : eccentricity,
    'y1': y1, 'x1': x1, 'Dx': Dx, 'y3': y3, 'yP2': yP2, 'fP3': fP3,
    'x2_derived' : x2_derived,
    'x3_derived' : x3_derived,
    'xP3_derived': xP3_derived,
    'yP3_derived': yP3_derived,
}

existing_rows = []
if os.path.exists(RESULTS_CSV):
    try:
        with open(RESULTS_CSV, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                existing_rows.append(row)
    except Exception as e:
        print('WARNING: Could not read existing results.csv - %s' % str(e))

updated = False
for i, row in enumerate(existing_rows):
    if str(row.get('SPECIMEN_NUM', '')).strip() == str(csv_num):
        existing_rows[i] = new_row
        updated = True
        break
if not updated:
    existing_rows.append(new_row)

try:
    existing_rows.sort(key=lambda r: int(r.get('SPECIMEN_NUM', 0)))
except Exception:
    pass

# Python 3: open with newline='' so csv module controls line endings.
try:
    f = open(RESULTS_CSV, 'w', newline='')
except TypeError:
    # Python 2 fallback (in case anyone still runs this under Abaqus 2017).
    f = open(RESULTS_CSV, 'wb')
try:
    writer = csv.DictWriter(f, fieldnames=HEADER, lineterminator='\n')
    writer.writeheader()
    writer.writerows(existing_rows)
finally:
    f.close()

print('')
print('============================================')
print('  results.csv updated: %s' % RESULTS_CSV)
print('  SPECIMEN_NUM    : %s' % csv_num)
print('  MAX_FI_ALL      : %s' % max_fi_all)
print('  MAX_FI_PINCHING : %s' % max_fi_pinching)
print('  MAX_FI_FOLDING  : %s' % max_fi_folding)
print('  MAX_SE          : %s' % max_se)
print('  Eccentricity    : %s' % eccentricity)
print('  DVs: y1=%s x1=%s Dx=%s y3=%s yP2=%s fP3=%s' % (y1, x1, Dx, y3, yP2, fP3))
print('  Derived: x2=%s x3=%s xP3=%s yP3=%s' % (x2_derived, x3_derived, xP3_derived, yP3_derived))
print('============================================')
