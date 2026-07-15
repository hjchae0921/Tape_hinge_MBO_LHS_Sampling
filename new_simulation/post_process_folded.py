# -*- coding: utf-8 -*-
# =============================================================
#  post_process_folded.py   (Abaqus 2017 Py2 / Abaqus 2024 Py3 compatible)
#
#  Same FI formula as post_process_csv_ver3.py, but the MBO objective is now
#  the FULLY-FOLDED failure index: MAX_FI_FOLDED = max-over-elements of FI at
#  the LAST frame of the folding step.  MAX_FI_ALL / PINCHING / FOLDING (max
#  over all frames) are still recorded for reference, plus MAX_SE (max ALLSE).
#
#  Run via wrapper that sets csv_num + BASE_DIR before exec().
# =============================================================
from __future__ import print_function
import os
import csv
import math

from abaqus import mdb, session
from abaqusConstants import *
from math import sqrt as _sqrt

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
print('  Post (FOLDED objective) specimen %d' % csv_num)
print('  ODB: %s' % odbPath)
print('============================================')

# ---- read design variables from cutout csv ----
eccentricity = ''
y1 = ''; x1 = ''; Dx = ''; y3 = ''; yP2 = ''; fP3 = ''
x2_derived = ''; x3_derived = ''; xP3_derived = ''; yP3_derived = ''
try:
    with open(csv_specimen, 'r') as f:
        rows = list(csv.reader(f))
    if len(rows) > 1:
        dr = rows[1]
        if len(dr) >= 3 and dr[2].strip():  eccentricity = float(dr[2].strip())
        if len(dr) >= 9:
            if dr[3].strip(): y1  = float(dr[3].strip())
            if dr[4].strip(): x1  = float(dr[4].strip())
            if dr[5].strip(): Dx  = float(dr[5].strip())
            if dr[6].strip(): y3  = float(dr[6].strip())
            if dr[7].strip(): yP2 = float(dr[7].strip())
            if dr[8].strip(): fP3 = float(dr[8].strip())
        if len(dr) >= 13:
            if dr[9].strip():  x2_derived  = float(dr[9].strip())
            if dr[10].strip(): x3_derived  = float(dr[10].strip())
            if dr[11].strip(): xP3_derived = float(dr[11].strip())
            if dr[12].strip(): yP3_derived = float(dr[12].strip())
except Exception as e:
    print('WARNING: failed to read cutout CSV - %s' % str(e))

# ---- open ODB ----
if odbPath in session.odbs.keys():
    session.odbs[odbPath].close()
odb = session.openOdb(name=odbPath, readOnly=True)
inst_keys = [k for k in odb.rootAssembly.instances.keys() if k.upper() != 'ASSEMBLY']
inst_obj = odb.rootAssembly.instances[inst_keys[0]]
print('Using instance: %s' % inst_keys[0])

# ---- failure-strength parameters (identical to production post) ----
sF1t = 139.47;  sF1c = 63.42
sF3  = 17.73;   sF4  = 5.07;   sF6 = 1.53
sJ1  = 1.0/sF1t - 1.0/sF1c
sK11 = 1.0/(sF1t*sF1c)
sK33 = 1.0/(sF3**2)
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


def frame_max_fi(frame):
    """Model-max FI at one frame (or None if SF/SM missing)."""
    available = list(frame.fieldOutputs.keys())
    if 'SF' not in available or 'SM' not in available:
        return None
    fValues = frame.fieldOutputs['SF'].getSubset(region=inst_obj).values
    mValues = frame.fieldOutputs['SM'].getSubset(region=inst_obj).values
    mMap = {}
    for mv in mValues:
        mMap[(int(mv.elementLabel), int(mv.integrationPoint))] = mv.data
    fmax = -1e99
    for fv in fValues:
        key = (int(fv.elementLabel), int(fv.integrationPoint))
        if key not in mMap:
            continue
        N1, N2, N3, N12, N13, N23 = fv.data
        M1, M2, M12 = mMap[key]
        Nx  = 0.5*(N1 + N2) + N12
        Ny  = 0.5*(N1 + N2) - N12
        Nxy = 0.5*(N2 - N1)
        Mx  = 0.5*(M1 + M2) + M12
        My  = 0.5*(M1 + M2) - M12
        M = Mx if abs(Mx) >= abs(My) else My
        fIndexIP = sJ1*(Nx+Ny) + sK11*(Nx**2 + Ny**2) + sK12*Nx*Ny + sK33*Nxy**2
        if fIndexIP >= 1.0:
            fIndexC = 1.0 + abs(M)/sF4
        else:
            dx = (sJ1 + sK12*Ny)**2 - 4*sK11*(sJ1*Ny + sK11*Ny**2 + sK33*Nxy**2 - 1)
            if dx >= 0:
                sd = _sqrt(dx)
                sFx = _select_root(Nx, (-(sJ1+sK12*Ny)+sd)/(2*sK11),
                                       (-(sJ1+sK12*Ny)-sd)/(2*sK11))
            else:
                sFx = None
            if sFx is None or abs(sFx) < 1e-9:
                fIndexCx = fIndexIP + abs(M)/sF4
            else:
                fIndexCx = abs(Nx/sFx) + abs(M)/sF4
            dy = (sJ1 + sK12*Nx)**2 - 4*sK11*(sJ1*Nx + sK11*Nx**2 + sK33*Nxy**2 - 1)
            if dy >= 0:
                sd = _sqrt(dy)
                sFy = _select_root(Ny, (-(sJ1+sK12*Nx)+sd)/(2*sK11),
                                       (-(sJ1+sK12*Nx)-sd)/(2*sK11))
            else:
                sFy = None
            if sFy is None or abs(sFy) < 1e-9:
                fIndexCy = fIndexIP + abs(M)/sF4
            else:
                fIndexCy = abs(Ny/sFy) + abs(M)/sF4
            fIndexC = fIndexCx if fIndexCx >= fIndexCy else fIndexCy
        if fIndexC > fmax:
            fmax = fIndexC
    return fmax if fmax > -1e99 else None


# ---- per-step max + fully-folded (last folding frame) ----
available_steps = list(odb.steps.keys())
stepNames = [s for s in ['pinching', 'folding'] if s in available_steps]

max_fi_all = -1e99; max_fi_pinching = -1e99; max_fi_folding = -1e99
max_fi_folded = ''      # objective: FI at last folding frame
for step_name in stepNames:
    step = odb.steps[step_name]
    step_max = -1e99
    n = len(step.frames)
    for fi_idx, frame in enumerate(step.frames):
        fm = frame_max_fi(frame)
        if fm is None:
            continue
        if fm > step_max:
            step_max = fm
        if step_name == 'folding' and fi_idx == n - 1:
            max_fi_folded = fm
    if step_name == 'pinching' and step_max > -1e99: max_fi_pinching = step_max
    if step_name == 'folding' and step_max > -1e99:  max_fi_folding = step_max
    if step_max > max_fi_all: max_fi_all = step_max
    print('  step %s max=%.6f' % (step_name, step_max if step_max > -1e99 else 0))

if max_fi_all      <= -1e99: max_fi_all      = ''
if max_fi_pinching <= -1e99: max_fi_pinching = ''
if max_fi_folding  <= -1e99: max_fi_folding  = ''
print('MAX_FI_FOLDED (objective) = %s' % max_fi_folded)

# ---- MAX_SE ----
max_se = ''
try:
    all_se = []
    for step_name in stepNames:
        step = odb.steps[step_name]
        for rk in step.historyRegions.keys():
            reg = step.historyRegions[rk]
            if 'ALLSE' in reg.historyOutputs:
                for tv, dv in reg.historyOutputs['ALLSE'].data:
                    all_se.append(dv)
    if all_se:
        max_se = max(all_se)
        print('MAX_SE = %s' % max_se)
    else:
        print('WARNING: ALLSE not found.')
except Exception as e:
    print('WARNING: MAX_SE failed - %s' % str(e))

session.odbs[odbPath].close()

# ---- write/update results.csv ----
HEADER = [
    'SPECIMEN_NUM',
    'MAX_FI_ALL', 'MAX_FI_PINCHING', 'MAX_FI_FOLDING', 'MAX_FI_FOLDED',
    'MAX_SE',
    'Eccentricity', 'y1', 'x1', 'Dx', 'y3', 'yP2', 'fP3',
    'x2_derived', 'x3_derived', 'xP3_derived', 'yP3_derived',
]
new_row = {
    'SPECIMEN_NUM': csv_num,
    'MAX_FI_ALL': max_fi_all, 'MAX_FI_PINCHING': max_fi_pinching,
    'MAX_FI_FOLDING': max_fi_folding, 'MAX_FI_FOLDED': max_fi_folded,
    'MAX_SE': max_se, 'Eccentricity': eccentricity,
    'y1': y1, 'x1': x1, 'Dx': Dx, 'y3': y3, 'yP2': yP2, 'fP3': fP3,
    'x2_derived': x2_derived, 'x3_derived': x3_derived,
    'xP3_derived': xP3_derived, 'yP3_derived': yP3_derived,
}

existing = []
if os.path.exists(RESULTS_CSV):
    try:
        with open(RESULTS_CSV, 'r') as f:
            existing = list(csv.DictReader(f))
    except Exception as e:
        print('WARNING: read results.csv - %s' % str(e))

updated = False
for i, row in enumerate(existing):
    if str(row.get('SPECIMEN_NUM', '')).strip() == str(csv_num):
        existing[i] = new_row; updated = True; break
if not updated:
    existing.append(new_row)
try:
    existing.sort(key=lambda r: int(r.get('SPECIMEN_NUM', 0)))
except Exception:
    pass

# Py3: newline='' ; Py2 fallback: 'wb'
try:
    f = open(RESULTS_CSV, 'w', newline='')
except TypeError:
    f = open(RESULTS_CSV, 'wb')
try:
    w = csv.DictWriter(f, fieldnames=HEADER, lineterminator='\n', extrasaction='ignore')
    w.writeheader(); w.writerows(existing)
finally:
    f.close()

print('results.csv updated: SPEC %s  FOLDED_FI=%s  SE=%s'
      % (csv_num, max_fi_folded, max_se))
