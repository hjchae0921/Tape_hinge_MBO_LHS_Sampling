#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  MODULE 4: STEP / INTERACTION / LOAD
#  - Amplitudes, Steps (pinching + folding)
#  - General Contact (ALL EXTERIOR)
#  - Boundary Conditions (shell-node BCs + rotation BCs)
#  - Loads (viscous pressure + concentrated pinching forces)
#  - Prerequisite: 01_part.py ~ 03_assembly.py
#  - Run: abaqus cae noGUI=04_step_interaction_load.py
# =============================================================
from abaqus import *
from abaqusConstants import *
from step import *
from interaction import *
from load import *
import regionToolset
import math

L_total = 220.0
D       = 38.0
R       = D / 2.0
PINCH_OFFSET = 5.0
FOLD_ANGLE   = 1.4835  # rad per end (~85 deg, total ~170 deg)

myModel = mdb.models['Model-1']
myPart  = myModel.parts['TapeSpring_3D']
a       = myModel.rootAssembly
inst    = a.instances['TapeSpring_3D-1']

# Recompute R_x, R_y from part geometry
# (conservative bounds for vertex selection)
R_x = 25.0
R_y = 25.0

# =============================================================
# [8] AMPLITUDES + STEPS
# =============================================================
# Pinching amplitude (TOTAL TIME, SMOOTH STEP)
#   0->1s: ramp up to 1.0,  1->4s: hold at 1.0
myModel.SmoothStepAmplitude(
    name='Pinching',
    timeSpan=TOTAL,
    data=((0.0, 0.0), (1.0, 1.0), (4.0, 1.0), (5.0, 0.0), (8.0, 0.0)))

# Rotation amplitude (TOTAL TIME, SMOOTH STEP)
#   0->1s: hold at 0.0,  1->4s: ramp up to 1.0
myModel.SmoothStepAmplitude(
    name='Rotation',
    timeSpan=TOTAL,
    data=((0.0, 0.0), (1.0, 0.0), (4.0, 1.0), (5.0, 1.0), (8.0, 0.0)))

# --- Step 1: Pinching (1 s, fixed mass scaling) ---
myModel.ExplicitDynamicsStep(
    name='pinching',
    previous='Initial',
    timePeriod=1.0,
    linearBulkViscosity=0.06,
    quadBulkViscosity=1.2,
    massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                  BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

# --- Step 2: Folding (3 s, variable mass scaling freq=1) ---
myModel.ExplicitDynamicsStep(
    name='folding',
    previous='pinching',
    timePeriod=3.0,
    linearBulkViscosity=0.06,
    quadBulkViscosity=1.2,
    massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                  BELOW_MIN, 1, 0, 0.0, 0.0, 1, None), ))

# --- Field Output ---
try:
    myModel.fieldOutputRequests['F-Output-1'].setValues(
        timeInterval=0.01,
        variables=('CF', 'RF', 'RM', 'RT', 'SE', 'SF', 'U', 'UR', 'UT'))
except KeyError:
    myModel.FieldOutputRequest('F-Output-1',
        createStepName='pinching', timeInterval=0.01,
        variables=('CF', 'RF', 'RM', 'RT', 'SE', 'SF', 'U', 'UR', 'UT'))

# --- History Output ---
try:
    myModel.historyOutputRequests['H-Output-1'].setValues(
        timeInterval=0.01,
        variables=('ALLAE', 'ALLCD', 'ALLCW', 'ALLDC', 'ALLDMD', 'ALLFD',
                   'ALLIE', 'ALLKE', 'ALLMW', 'ALLPD', 'ALLPW', 'ALLSE',
                   'ALLVD', 'ALLWK'))
except KeyError:
    myModel.HistoryOutputRequest('H-Output-1',
        createStepName='pinching', timeInterval=0.01,
        variables=('ALLAE', 'ALLCD', 'ALLCW', 'ALLDC', 'ALLDMD', 'ALLFD',
                   'ALLIE', 'ALLKE', 'ALLMW', 'ALLPD', 'ALLPW', 'ALLSE',
                   'ALLVD', 'ALLWK'))

print('===== Steps created =====')

# =============================================================
# [9] INTERACTION MODULE
# =============================================================
myModel.ContactProperty('IntProp-1')
myModel.interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
myModel.interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON,
    constraintEnforcementMethod=DEFAULT,
    pressureOverclosure=HARD)

myModel.ContactExp(createStepName='Initial', name='Int-1')
myModel.interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
myModel.interactions['Int-1'].contactPropertyAssignments.appendInStep(
    assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

print('===== Interaction created =====')

# =============================================================
# [10] BOUNDARY CONDITIONS + LOADS
# =============================================================
tol_bc = 0.5

# --- Initial BCs on shell nodes (prevent rigid-body motion) ---
# BC-1: Y-extreme vertices at Z=L_total -> u1=0
v_bc1 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc, yMin=-R_y - tol_bc, zMin=L_total - tol_bc,
    xMax= tol_bc, yMax= R_y + tol_bc, zMax=L_total + tol_bc)
a.Set(name='BC1-nodes', vertices=v_bc1)
print('BC-1: %d vertices (Y-extreme at Z=L)' % len(v_bc1))

myModel.DisplacementBC(
    name='BC-1', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC1-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# BC-2: X-extreme vertices at Z=L_total -> u2=0, u3=0
v_bc2 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc, zMin=L_total - tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc, zMax=L_total + tol_bc)
a.Set(name='BC2-nodes', vertices=v_bc2)
print('BC-2: %d vertices (X-extreme at Z=L)' % len(v_bc2))

myModel.DisplacementBC(
    name='BC-2', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC2-nodes'],
    u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# BC-3: Y-extreme vertices at Z=0 -> u1=0
v_bc3 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc, yMin=-R_y - tol_bc, zMin=-tol_bc,
    xMax= tol_bc, yMax= R_y + tol_bc, zMax= tol_bc)
a.Set(name='BC3-nodes', vertices=v_bc3)
print('BC-3: %d vertices (Y-extreme at Z=0)' % len(v_bc3))

myModel.DisplacementBC(
    name='BC-3', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC3-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# BC-4: X-extreme vertices at Z=0 -> u2=0
v_bc4 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc, zMin=-tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc, zMax= tol_bc)
a.Set(name='BC4-nodes', vertices=v_bc4)
print('BC-4: %d vertices (X-extreme at Z=0)' % len(v_bc4))

myModel.DisplacementBC(
    name='BC-4', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC4-nodes'],
    u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# --- Viscous pressure on all faces (created in pinching, propagates) ---
allFaces = inst.faces[:]
myModel.Pressure(
    amplitude=UNSET,
    createStepName='pinching',
    distributionType=VISCOUS,
    field='',
    magnitude=1.4623e-06,
    name='Load-1',
    region=Region(side1Faces=allFaces))

print('===== Viscous pressure applied =====')

# --- Concentrated pinching forces ---
# 9 top vertices (Y > 0) at X={-5,0,5} x Z={120,125,130}
# 9 bottom vertices (Y < 0) at X={-5,0,5} x Z={120,125,130}
z_mid = L_total / 2.0

# Compute the minimum |Y| on the ellipse at X = PINCH_OFFSET
# For R_x, R_y we need exact values. Recompute from part name/geometry.
# Safe approach: select top vertices by large positive Y threshold
# Top of cylinder: Y ranges from ~R_y*sqrt(1-(5/R_x)^2) to R_y
# Use generous lower bound for Y
y_threshold_top = 10.0   # well above cutout, below cylinder crown

v_pinch_top = inst.vertices.getByBoundingBox(
    xMin=-PINCH_OFFSET - tol_bc,
    yMin= y_threshold_top,
    zMin= z_mid - PINCH_OFFSET - tol_bc,
    xMax= PINCH_OFFSET + tol_bc,
    yMax= R_y + tol_bc,
    zMax= z_mid + PINCH_OFFSET + tol_bc)

v_pinch_bot = inst.vertices.getByBoundingBox(
    xMin=-PINCH_OFFSET - tol_bc,
    yMin=-R_y - tol_bc,
    zMin= z_mid - PINCH_OFFSET - tol_bc,
    xMax= PINCH_OFFSET + tol_bc,
    yMax=-y_threshold_top,
    zMax= z_mid + PINCH_OFFSET + tol_bc)

a.Set(name='pinch-top', vertices=v_pinch_top)
a.Set(name='pinch-bot', vertices=v_pinch_bot)

print('Pinching vertices: top=%d, bottom=%d (expected 9 each)' %
      (len(v_pinch_top), len(v_pinch_bot)))

# Top pinching: push downward (CF2 = -3.0)
myModel.ConcentratedForce(
    name='pinchload-top',
    createStepName='pinching',
    region=a.sets['pinch-top'],
    cf1=0.0, cf2=-3.0, cf3=0.0,
    amplitude='Pinching')

# Bottom pinching: push upward (CF2 = +3.0)
myModel.ConcentratedForce(
    name='pinchload-bot',
    createStepName='pinching',
    region=a.sets['pinch-bot'],
    cf1=0.0, cf2=3.0, cf3=0.0,
    amplitude='Pinching')

print('===== Pinching loads applied =====')

# --- Folding rotation BCs (created in folding step) ---
myModel.DisplacementBC(
    name='BC-5', createStepName='folding',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['RP-A'],
    u1=UNSET, u2=UNSET, u3=UNSET,
    ur1=FOLD_ANGLE,
    ur2=UNSET, ur3=UNSET,
    amplitude='Rotation')

myModel.DisplacementBC(
    name='BC-6', createStepName='folding',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['RP-B'],
    u1=UNSET, u2=UNSET, u3=UNSET,
    ur1=-FOLD_ANGLE,
    ur2=UNSET, ur3=UNSET,
    amplitude='Rotation')

print('===== Folding rotation BCs applied =====')
print('===== Step/Interaction/Load module completed =====')
