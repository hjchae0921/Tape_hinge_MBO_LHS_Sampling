#! /user/bin/python
# -*- coding: UTF-8 -*-
# -*- coding: mbcs -*-
# get_FI_graph_pinching.py  [FIXED VERSION]
# Extract max FI1 (in-plane), FI2 (moment), FI3 (interaction) vs time
# from the currently open ODB (pinching 1s + folding 3s).
# X-axis: Total Time (s),  Y-axis: Failure Index
# Save plot as PNG + CSV data.
# Usage: Run via File -> Run Script in Abaqus CAE with an ODB already open.
#
# === FIX SUMMARY (vs original) ===
# v1 fix: branch on sign of Nx/Ny (not Ny/Nx); compute both roots and
#         select same-sign root; fall back to FI1 on degenerate case.
# v2 fix: preserve bending term abs(M)/sF4 even when FI1 >= 1.
#         Without this, FI3 would drop discontinuously from ~1.17 to
#         ~1.0 each time FI1 crosses 1 during folding, which is
#         non-physical (bending keeps acting).

from abaqus import *
from abaqusConstants import *
from odbAccess import *
import visualization
import csv
import os
import math

# ============================================================
# Failure strength parameters (must match post_process_FI.py)
# ============================================================
sF1t = 139.47
sF1c = 63.42
sF3  = 17.73
sF4  = 5.07
sF6  = 1.53

sJ1  = 1.0/sF1t - 1.0/sF1c
sK11 = 1.0/(sF1t * sF1c)
sK33 = 1.0/(sF3**2)
sK44 = 1.0/(sF4**2)
sK66 = 1.0/(sF6**2)
sK12 = -sK11 / 2.0


# ============================================================
# Helper: pick the physically correct root of Fx (or Fy)
# ============================================================
def _select_root(N_self, root_plus, root_minus):
    """
    Select the root with the same sign as N_self.
    If both roots have the same sign as N_self, pick the one
    closer to N_self (the actual failure boundary on this side).
    If neither matches in sign, return None (degenerate case).
    """
    # Same-sign filter
    if N_self >= 0:
        candidates = [r for r in (root_plus, root_minus) if r > 0]
    else:
        candidates = [r for r in (root_plus, root_minus) if r < 0]

    if not candidates:
        return None  # No root on the same side -> degenerate

    if len(candidates) == 1:
        return candidates[0]

    # Both roots on the same side: pick the one with larger magnitude
    # (the outer boundary, which is the actual ultimate strength)
    return candidates[0] if abs(candidates[0]) >= abs(candidates[1]) else candidates[1]


# ============================================================
# 0. Detect currently open ODB
# ============================================================
if len(session.odbs) == 0:
    raise RuntimeError('No ODB is open. Please open an ODB first.')

odb_name = session.odbs.keys()[0]
odb = session.odbs[odb_name]
print('Using ODB: %s' % odb_name)

# Get instance (skip ASSEMBLY)
inst_keys = [k for k in odb.rootAssembly.instances.keys()
             if k.upper() != 'ASSEMBLY']
if len(inst_keys) == 0:
    raise KeyError('No part instance found in ODB.')
inst_obj = odb.rootAssembly.instances[inst_keys[0]]
print('Using instance: %s' % inst_keys[0])

# ============================================================
# 1. Select steps: pinching (1s) + folding (3s)
# ============================================================
available_steps = list(odb.steps.keys())
print('Available steps in ODB: %s' % str(available_steps))

target_steps = ['pinching', 'folding']
stepNames = [s for s in target_steps if s in available_steps]
if not stepNames:
    stepNames = [s for s in available_steps
                 if len(odb.steps[s].frames) > 1]
print('Processing steps: %s' % str(stepNames))

for sn in stepNames:
    s = odb.steps[sn]
    print('  Step "%s": %d frames, timePeriod=%.2f' % (sn, len(s.frames), s.timePeriod))

# ============================================================
# 2. Process all frames across steps
# ============================================================
results = []
_sqrt = math.sqrt
timeOffset = 0.0

for sName in stepNames:
    step = odb.steps[sName]
    nFrames = len(step.frames)
    print('--- Processing step: %s (%d frames) ---' % (sName, nFrames))

    for frame_idx in range(nFrames):
        try:
            frame = step.frames[frame_idx]
            available = list(frame.fieldOutputs.keys())
        except Exception as e:
            print('  Frame %d: skipped (corrupt) - %s' % (frame_idx, str(e)))
            continue

        if 'SF' not in available or 'SM' not in available:
            continue

        totalTime = timeOffset + frame.frameValue

        try:
            fValues = frame.fieldOutputs['SF'].getSubset(region=inst_obj).values
            mValues = frame.fieldOutputs['SM'].getSubset(region=inst_obj).values
        except Exception as e:
            print('  Frame %d: skipped (field read error) - %s' % (frame_idx, str(e)))
            continue

        # Build SM lookup by (elementLabel, integrationPoint)
        mMap = {}
        for mv in mValues:
            mMap[(int(mv.elementLabel), int(mv.integrationPoint))] = mv.data

        maxFI1 = 0.0
        maxFI2 = 0.0
        maxFI3 = 0.0

        for fv in fValues:
            key = (int(fv.elementLabel), int(fv.integrationPoint))
            if key not in mMap:
                continue

            N1, N2, N3, N12, N13, N23 = fv.data
            M1, M2, M12 = mMap[key]

            Nx  =  0.5*(N1 + N2) + N12
            Ny  =  0.5*(N1 + N2) - N12
            Nxy =  0.5*(N2 - N1)
            Mx  =  0.5*(M1 + M2) + M12
            My  =  0.5*(M1 + M2) - M12
            Mxy =  0.5*(M2 - M1)

            M = Mx if abs(Mx) >= abs(My) else My

            # FI1: In-plane
            fIndexIP = (sJ1*(Nx + Ny)
                        + sK11*(Nx**2 + Ny**2)
                        + sK12*Nx*Ny
                        + sK33*Nxy**2)

            # FI2: Moment
            fIndexM = sK44*M**2 + sK66*Mxy**2

            # ============================================================
            # FI3: Interaction (FIXED v2)
            # FIX v2: when FI1>=1, preserve bending term to avoid the
            # spurious drop in FI3 (it would otherwise step down from
            # ~1.17 to ~1.0 the moment FI1 crosses 1, even though
            # physically bending is still acting).
            # ============================================================
            if fIndexIP >= 1.0:
                # Outside FI1 ellipse: clamp in-plane ratio to 1,
                # but keep the bending contribution.
                fIndexC = 1.0 + abs(M)/sF4
            else:
                # --- Fx: Nx-direction failure boundary on FI1=1 ellipse
                dx = ((sJ1 + sK12*Ny)**2
                      - 4*sK11*(sJ1*Ny + sK11*Ny**2 + sK33*Nxy**2 - 1))

                if dx >= 0:
                    sqrt_dx = _sqrt(dx)
                    Fx_plus  = (-(sJ1 + sK12*Ny) + sqrt_dx) / (2*sK11)
                    Fx_minus = (-(sJ1 + sK12*Ny) - sqrt_dx) / (2*sK11)
                    sFx = _select_root(Nx, Fx_plus, Fx_minus)
                else:
                    sFx = None

                if sFx is None or abs(sFx) < 1e-9:
                    # Degenerate fallback: include bending term
                    fIndexCx = fIndexIP + abs(M)/sF4
                else:
                    fIndexCx = abs(Nx/sFx) + abs(M)/sF4

                # --- Fy: Ny-direction failure boundary on FI1=1 ellipse
                dy = ((sJ1 + sK12*Nx)**2
                      - 4*sK11*(sJ1*Nx + sK11*Nx**2 + sK33*Nxy**2 - 1))

                if dy >= 0:
                    sqrt_dy = _sqrt(dy)
                    Fy_plus  = (-(sJ1 + sK12*Nx) + sqrt_dy) / (2*sK11)
                    Fy_minus = (-(sJ1 + sK12*Nx) - sqrt_dy) / (2*sK11)
                    sFy = _select_root(Ny, Fy_plus, Fy_minus)
                else:
                    sFy = None

                if sFy is None or abs(sFy) < 1e-9:
                    fIndexCy = fIndexIP + abs(M)/sF4
                else:
                    fIndexCy = abs(Ny/sFy) + abs(M)/sF4

                fIndexC = fIndexCx if fIndexCx >= fIndexCy else fIndexCy

            if fIndexIP > maxFI1:
                maxFI1 = fIndexIP
            if fIndexM > maxFI2:
                maxFI2 = fIndexM
            if fIndexC > maxFI3:
                maxFI3 = fIndexC

        results.append((totalTime, sName, maxFI1, maxFI2, maxFI3))

        if frame_idx % 50 == 0:
            print('  Frame %d/%d: t=%.3f s, FI1=%.4f, FI2=%.4f, FI3=%.4f'
                  % (frame_idx, nFrames, totalTime, maxFI1, maxFI2, maxFI3))

    timeOffset += step.timePeriod

print('Total frames processed: %d' % len(results))

if len(results) == 0:
    raise RuntimeError('No valid frames with SF/SM data found. Check ODB.')

# ============================================================
# 3. Save CSV
# ============================================================
out_dir = os.path.dirname(odb_name)
csv_path = os.path.join(out_dir, 'FI_graph_data.csv')
with open(csv_path, 'w') as f:
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow(['TotalTime_s', 'Step', 'max_FI1', 'max_FI2', 'max_FI3'])
    for row in results:
        writer.writerow(['%.6f' % row[0], row[1],
                         '%.6f' % row[2], '%.6f' % row[3], '%.6f' % row[4]])
print('CSV saved: %s' % csv_path)

# ============================================================
# 4. Plot using Abaqus XY Plot and save as image
# ============================================================
times    = [r[0] for r in results]
fi1_vals = [r[2] for r in results]
fi2_vals = [r[3] for r in results]
fi3_vals = [r[4] for r in results]

# Clean up old XYData if exists
for nm in ['max_FI1', 'max_FI2', 'max_FI3']:
    if nm in session.xyDataObjects:
        del session.xyDataObjects[nm]

xy_fi1 = session.XYData(
    name='max_FI1',
    data=tuple(zip(times, fi1_vals)),
    xValuesLabel='Time (s)',
    yValuesLabel='Failure Index',
    legendLabel='max FI1 (In-plane)')

xy_fi2 = session.XYData(
    name='max_FI2',
    data=tuple(zip(times, fi2_vals)),
    xValuesLabel='Time (s)',
    yValuesLabel='Failure Index',
    legendLabel='max FI2 (Moment)')

xy_fi3 = session.XYData(
    name='max_FI3',
    data=tuple(zip(times, fi3_vals)),
    xValuesLabel='Time (s)',
    yValuesLabel='Failure Index',
    legendLabel='max FI3 (Interaction)')

c1 = session.Curve(xyData=xy_fi1)
c2 = session.Curve(xyData=xy_fi2)
c3 = session.Curve(xyData=xy_fi3)

plotName = 'FI_vs_Time'
if plotName in session.xyPlots:
    del session.xyPlots[plotName]
xyPlot = session.XYPlot(plotName)
chartName = xyPlot.charts.keys()[0]
chart = xyPlot.charts[chartName]
chart.setValues(curvesToPlot=(c1, c2, c3))

chart.axes1[0].axisData.setValues(title='Time (s)')
chart.axes2[0].axisData.setValues(title='Failure Index')
chart.legend.setValues(show=True)

vp = session.viewports['Viewport: 1']
vp.setValues(displayedObject=xyPlot)

# Save as PNG
img_path = os.path.join(out_dir, 'FI_graph')
try:
    session.printToFile(
        fileName=img_path, format=PNG,
        canvasObjects=(vp,))
    print('Graph saved: %s.png' % img_path)
except Exception as e:
    print('WARNING: Could not save graph image - %s' % str(e))
    print('Data is available in CSV: %s' % csv_path)

print('===== Done =====')