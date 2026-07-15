#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  TapeSpring Hinge — Pinching + Folding (all-in-one)
#  Part -> Property -> Assembly -> Step/Interaction/Load -> Mesh -> Job
#
#  Output directory: E:\Tape_BO\pinching\sim_<csv_num>\
#    All job files (INP, ODB, STA, MSG, DAT ...) are written there.
#
#  Run:  abaqus cae noGUI=tube_hinge_pinching.py
#  Compatible with Abaqus 2017 (Python 2.7) and 2024 (Python 3).
# =============================================================
from __future__ import print_function
from abaqus import *
from abaqusConstants import *
import __main__
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import regionToolset
import sys, os
import math
import csv

# =============================================================
# [0] Parameters
# =============================================================
L_total = 220.0          # Hinge length (Z-axis), mm
D       = 38.0           # Cylinder diameter, mm
R       = D / 2.0        # Radius = 19.0 mm
try:
    csv_num
except NameError:
    csv_num = 10

COUPLING_STRIP = 5.0     # mm from each end for coupling surface
PINCH_OFFSET   = 5.0     # mm offset from center for pinching grid
PINCH_FORCE    = 3.0     # N, concentrated force magnitude per node
FOLD_ANGLE     = 1.4835  # rad per end (~85 deg, total relative ~170 deg)

# Repo-relative paths: BASE_DIR is runner-injected; else derived from this script's location.
try:
    BASE_DIR
except NameError:
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        BASE_DIR = os.getcwd()

SPECIMEN_DIR = os.path.join(BASE_DIR, 'specimen')
SIM_DIR      = os.path.join(BASE_DIR, 'sim_%d' % csv_num)
csv_file     = os.path.join(SPECIMEN_DIR, 'cutout_points_%d.csv' % csv_num)

myModel = mdb.models['Model-1']

# =============================================================
# [1] Read 2D cutout coordinates and eccentricity from CSV
# =============================================================
points_2d = []
eccentricity = 0.0

with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        if i == 0 and row[0].strip().upper() == 'X':
            continue
        if len(row) < 2:
            continue
        a, b = row[0].strip(), row[1].strip()
        if a and b:
            points_2d.append((float(a), float(b)))
        if len(row) >= 3 and row[2].strip() and eccentricity == 0.0:
            try:
                eccentricity = float(row[2].strip())
            except ValueError:
                pass

print('Loaded %d points from CSV. Eccentricity: %f' % (len(points_2d), eccentricity))

# =============================================================
# [2] PART MODULE — Base Elliptical Shell
# =============================================================
k = math.sqrt(1.0 - eccentricity**2)
denom = 3.0 * (1.0 + k) - math.sqrt((3.0 + k) * (1.0 + 3.0 * k))
R_x = 2.0 * R / denom
R_y = R_x * k

s_base = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
s_base.EllipseByCenterPerimeter(center=(0.0, 0.0),
                                axisPoint1=(R_x, 0.0),
                                axisPoint2=(0.0, R_y))

myPart = myModel.Part(name='TapeSpring_3D',
                      dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)
myPart.BaseShellExtrude(sketch=s_base, depth=L_total)
del myModel.sketches['__profile__']
print('Base elliptical shell (Rx=%.4f, Ry=%.4f, L=%.1f)' % (R_x, R_y, L_total))

# =============================================================
# [3] 2D -> 3D projection
# =============================================================
# Cutout point order from bezier_c2_cutout_sampling/mirror_to_full():
#   TL(100) -> TR(100) -> BR[1:](99) -> BL[1:](99) -> close(1) = 399 pts
# Starts at left-axis crossing (-x3, 0), traverses CCW, closes back to start.
# Find the two unique Y=0 axis crossings, then split into top/bot halves (both LTR).

# Drop closing duplicate (last point == first point)
if len(points_2d) >= 2 and points_2d[0] == points_2d[-1]:
    points_2d = points_2d[:-1]

zero_idx = [i for i in range(len(points_2d)) if abs(points_2d[i][1]) < 1e-6]
if len(zero_idx) >= 2:
    i_left  = zero_idx[0]   # left  axis crossing (-x3, 0)
    i_right = zero_idx[1]   # right axis crossing (+x3, 0)
    # Top half LTR: i_left -> top -> i_right
    top_2d = points_2d[i_left : i_right + 1]
    # Bot half CCW (RTL): i_right -> bottom -> wrap back to i_left
    bot_2d_ccw = points_2d[i_right:] + [points_2d[i_left]]
    # Reverse to LTR to match the script's downstream convention
    bot_2d = list(reversed(bot_2d_ccw))
else:
    # Fallback (shouldn't happen for well-formed Bezier cutouts)
    top_2d = [p for p in points_2d if p[1] >= 0]
    bot_2d = [p for p in points_2d if p[1] <= 0]

top_pts = []
for x2d, y2d in top_2d:
    u = y2d
    v = x2d + (L_total / 2.0)
    top_pts.append((u, v))

bot_pts = []
for x2d, y2d in bot_2d:
    u = y2d
    v = x2d + (L_total / 2.0)
    bot_pts.append((u, v))

bot_pts.reverse()


def clean_pts(pts):
    clean = [pts[0]]
    for i in range(1, len(pts)):
        du = pts[i][0] - clean[-1][0]
        dv = pts[i][1] - clean[-1][1]
        if (du*du + dv*dv) > 1e-6:
            clean.append(pts[i])
    return clean


top_pts = clean_pts(top_pts)
bot_pts = clean_pts(bot_pts)
print('Cutout: %d top pts, %d bot pts' % (len(top_pts), len(bot_pts)))

# =============================================================
# [4] Datum Plane (YZ) + sketch + CutExtrude (both +-X)
# =============================================================
def draw_half_curve(sketch, sk_pts, y_tol=0.01):
    n = len(sk_pts)
    if n < 2:
        return
    y_vals = [p[1] for p in sk_pts]
    best_start, best_end, best_len = 0, 0, 0
    i = 0
    while i < n:
        y_ref = y_vals[i]
        j = i + 1
        while j < n and abs(y_vals[j] - y_ref) < y_tol:
            j += 1
        if (j - i) > best_len:
            best_len = j - i
            best_start = i
            best_end = j - 1
        i = j
    if best_len < 5:
        sketch.Spline(points=tuple(sk_pts))
        print('  -> single Spline (%d pts)' % n)
        return
    left_curve  = sk_pts[:best_start + 1]
    right_curve = sk_pts[best_end:]
    if len(left_curve) >= 2:
        sketch.Spline(points=tuple(left_curve))
    sketch.Line(point1=sk_pts[best_start], point2=sk_pts[best_end])
    if len(right_curve) >= 2:
        sketch.Spline(points=tuple(right_curve))
    print('  -> L-curve(%d) + Line + R-curve(%d)' %
          (len(left_curve), len(right_curve)))


dp_id = myPart.DatumPlaneByPrincipalPlane(
    principalPlane=YZPLANE, offset=0.0).id
da_id = myPart.DatumAxisByPrincipalAxis(
    principalAxis=ZAXIS).id

t = myPart.MakeSketchTransform(
    sketchPlane=myPart.datums[dp_id],
    sketchUpEdge=myPart.datums[da_id],
    sketchPlaneSide=SIDE1,
    sketchOrientation=TOP,
    origin=(0.0, 0.0, L_total / 2.0))

s_cut = myModel.ConstrainedSketch(
    name='__cutout__', sheetSize=500.0, transform=t)
myPart.projectReferencesOntoSketch(sketch=s_cut, filter=COPLANAR_EDGES)

top_sk = [(v - L_total / 2.0, u) for (u, v) in top_pts]
bot_sk = [(v - L_total / 2.0, u) for (u, v) in bot_pts]

if len(top_sk) > 0 and len(bot_sk) > 0:
    bot_sk[0]  = top_sk[-1]
    bot_sk[-1] = top_sk[0]

print('Drawing top half:')
draw_half_curve(s_cut, top_sk)
print('Drawing bot half:')
draw_half_curve(s_cut, bot_sk)

myPart.CutExtrude(
    sketchPlane=myPart.datums[dp_id],
    sketchUpEdge=myPart.datums[da_id],
    sketchPlaneSide=SIDE1,
    sketchOrientation=TOP,
    sketch=s_cut,
    flipExtrudeDirection=ON,
    depth=R_x + 5.0)
del myModel.sketches['__cutout__']
print('+X cutout applied.')

s_cut2 = myModel.ConstrainedSketch(
    name='__cutout2__', sheetSize=500.0, transform=t)
myPart.projectReferencesOntoSketch(sketch=s_cut2, filter=COPLANAR_EDGES)

print('Drawing top half (cut2):')
draw_half_curve(s_cut2, top_sk)
print('Drawing bot half (cut2):')
draw_half_curve(s_cut2, bot_sk)

myPart.CutExtrude(
    sketchPlane=myPart.datums[dp_id],
    sketchUpEdge=myPart.datums[da_id],
    sketchPlaneSide=SIDE1,
    sketchOrientation=TOP,
    sketch=s_cut2,
    flipExtrudeDirection=OFF,
    depth=R_x + 5.0)
del myModel.sketches['__cutout2__']
print('-X cutout applied.')

# =============================================================
# [5] Partitions
# =============================================================
# --- Basic: XY(Z=0), YZ(X=0), XZ(Y=0) ---
dp_xy = myPart.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
dp_yz = myPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
dp_xz = myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)

for dp_obj, name in [(dp_xy, 'XY(Z=0)'),
                      (dp_yz, 'YZ(X=0)'),
                      (dp_xz, 'XZ(Y=0)')]:
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp_obj.id],
            faces=myPart.faces[:])
        print('Partitioned at %s' % name)
    except Exception as e:
        print('%s partition skipped: %s' % (name, str(e)))

# --- Cutout bounding rectangle ---
all_cut_pts = top_pts + bot_pts
z_min_cut = min(v for (u, v) in all_cut_pts)
z_max_cut = max(v for (u, v) in all_cut_pts)
u_max_cut = max(abs(u) for (u, v) in all_cut_pts)

rect_spacing = 1.0
z_min_rect = z_min_cut - rect_spacing
z_max_rect = z_max_cut + rect_spacing
y_min_rect = -(u_max_cut + rect_spacing)
y_max_rect = u_max_cut + rect_spacing

print('Cutout bbox: Z=[%.1f, %.1f], Y=[%.1f, %.1f]' %
      (z_min_cut, z_max_cut, -u_max_cut, u_max_cut))

for offset_val, plane_type, label in [
    (z_min_rect, XYPLANE, 'Z=%.2f (rect bottom)'),
    (z_max_rect, XYPLANE, 'Z=%.2f (rect top)'),
    (y_min_rect, XZPLANE, 'Y=%.2f (rect left)'),
    (y_max_rect, XZPLANE, 'Y=%.2f (rect right)'),
]:
    dp = myPart.DatumPlaneByPrincipalPlane(
        principalPlane=plane_type, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id],
            faces=myPart.faces[:])
        print('Partitioned at %s' % (label % offset_val))
    except Exception as e:
        print('Partition at %s skipped: %s' % (label % offset_val, str(e)))

# --- Coupling strip: Z=5, Z=245 ---
for offset_val, label in [
    (COUPLING_STRIP,           'Z=%.1f (coupling bottom strip)'),
    (L_total - COUPLING_STRIP, 'Z=%.1f (coupling top strip)'),
]:
    dp = myPart.DatumPlaneByPrincipalPlane(
        principalPlane=XYPLANE, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id],
            faces=myPart.faces[:])
        print('Partitioned at %s' % (label % offset_val))
    except Exception as e:
        print('Partition at %s skipped: %s' % (label % offset_val, str(e)))

# --- Pinching grid: Z=120/125/130, X=+-5 ---
z_mid = L_total / 2.0
for offset_val, plane_type, label in [
    (z_mid - PINCH_OFFSET, XYPLANE, 'Z=%.1f (pinch bottom)'),
    (z_mid,                XYPLANE, 'Z=%.1f (pinch center)'),
    (z_mid + PINCH_OFFSET, XYPLANE, 'Z=%.1f (pinch top)'),
    (PINCH_OFFSET,         YZPLANE, 'X=%.1f (pinch +X)'),
    (-PINCH_OFFSET,        YZPLANE, 'X=%.1f (pinch -X)'),
]:
    dp = myPart.DatumPlaneByPrincipalPlane(
        principalPlane=plane_type, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id],
            faces=myPart.faces[:])
        print('Partitioned at %s' % (label % offset_val))
    except Exception as e:
        print('Partition at %s skipped: %s' % (label % offset_val, str(e)))

print('===== Part [TapeSpring_3D] created =====')

# =============================================================
# [6] PROPERTY MODULE
# =============================================================
myModel.GeneralStiffnessSection(
    applyThermalStress=0,
    density=3.18e-10,
    name='Section-1',
    poissonDefinition=DEFAULT,
    referenceTemperature=None,
    stiffnessMatrix=(
        7714.0, 6380.0, 7714.0, 0.0, 0.0, 5962.0,
        0.0, 0.0, 0.0, 23.6,
        0.0, 0.0, 0.0, 19.1, 23.6,
        0.0, 0.0, 0.0, 0.0, 0.0, 19.9),
    useDensity=ON)

myPart.SectionAssignment(
    offset=0.0,
    offsetField='',
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=myPart.faces[:]),
    sectionName='Section-1',
    thicknessAssignment=FROM_SECTION)

datumCsys = myPart.DatumCsysByThreePoints(
    coordSysType=CYLINDRICAL,
    name='Datum csys-cyl',
    origin=(0.0, 0.0, 0.0),
    point1=(0.0, 0.0, 1.0),
    point2=(1.0, 0.0, 0.0))

myPart.MaterialOrientation(
    additionalRotationField='',
    additionalRotationType=ROTATION_NONE,
    angle=0.0,
    axis=AXIS_3,
    fieldName='',
    localCsys=myPart.datums[datumCsys.id],
    orientationType=SYSTEM,
    region=Region(faces=myPart.faces[:]))

print('===== Property assigned =====')

# =============================================================
# [7] ASSEMBLY MODULE
# =============================================================
a = myModel.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(dependent=ON, name='TapeSpring_3D-1',
           part=myModel.parts['TapeSpring_3D'])

inst = a.instances['TapeSpring_3D-1']

# --- Reference Points ---
rpA_feat = a.ReferencePoint(point=(0.0, 0.0, L_total))
rpB_feat = a.ReferencePoint(point=(0.0, 0.0, 0.0))
rpC_feat = a.ReferencePoint(point=(0.0, 0.0, L_total / 2.0))

a.features.changeKey(fromName=rpA_feat.name, toName='RP-A')
a.features.changeKey(fromName=rpB_feat.name, toName='RP-B')
a.features.changeKey(fromName=rpC_feat.name, toName='RP-C')

rpA = a.referencePoints[rpA_feat.id]
rpB = a.referencePoints[rpB_feat.id]
rpC = a.referencePoints[rpC_feat.id]

a.Set(name='RP-A', referencePoints=(rpA,))
a.Set(name='RP-B', referencePoints=(rpB,))
a.Set(name='RP-C', referencePoints=(rpC,))

# --- Surface-based Kinematic Coupling ---
tol = 0.1

facesA = inst.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=-R_y - tol, zMin=L_total - COUPLING_STRIP - tol,
    xMax= R_x + tol, yMax= R_y + tol, zMax=L_total + tol)

facesB = inst.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=-R_y - tol, zMin=-tol,
    xMax= R_x + tol, yMax= R_y + tol, zMax=COUPLING_STRIP + tol)

print('Coupling faces: top=%d, bottom=%d' % (len(facesA), len(facesB)))

a.Surface(name='CoupSurf-Top', side1Faces=facesA)
a.Surface(name='CoupSurf-Bot', side1Faces=facesB)

myModel.Coupling(
    controlPoint=Region(referencePoints=(rpA,)),
    couplingType=KINEMATIC,
    influenceRadius=WHOLE_SURFACE,
    localCsys=None,
    name='Constraint-1',
    surface=a.surfaces['CoupSurf-Top'],
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

myModel.Coupling(
    controlPoint=Region(referencePoints=(rpB,)),
    couplingType=KINEMATIC,
    influenceRadius=WHOLE_SURFACE,
    localCsys=None,
    name='Constraint-2',
    surface=a.surfaces['CoupSurf-Bot'],
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

# --- Point Inertia on RP-C ---
a.engineeringFeatures.PointMassInertia(
    name='Inertia-1',
    region=a.sets['RP-C'],
    mass=1e-06,
    i11=1e-06, i22=1e-06, i33=1e-06)

print('===== Assembly + Coupling + Inertia done =====')

# =============================================================
# [8] STEP MODULE  (pinching 1 s + folding 3 s)
# =============================================================
# Optional folding-only mode (validation study): skip the pinching step so
# folding is applied directly. Injectable; defaults to False so normal
# production/MBO runs are byte-for-byte unchanged.
try:
    SKIP_PINCHING
except NameError:
    SKIP_PINCHING = False
FIRST_STEP = 'folding' if SKIP_PINCHING else 'pinching'

myModel.SmoothStepAmplitude(
    name='Pinching',
    timeSpan=TOTAL,
    data=((0.0, 0.0), (1.0, 1.0), (4.0, 1.0), (5.0, 0.0), (8.0, 0.0)))

if SKIP_PINCHING:
    # rotation ramps 0->1 over the single folding step (0..3 s)
    myModel.SmoothStepAmplitude(
        name='Rotation', timeSpan=TOTAL,
        data=((0.0, 0.0), (3.0, 1.0)))
else:
    myModel.SmoothStepAmplitude(
        name='Rotation', timeSpan=TOTAL,
        data=((0.0, 0.0), (1.0, 0.0), (4.0, 1.0), (5.0, 1.0), (8.0, 0.0)))

if not SKIP_PINCHING:
    myModel.ExplicitDynamicsStep(
        name='pinching',
        previous='Initial',
        timePeriod=1.0,
        linearBulkViscosity=0.06,
        quadBulkViscosity=1.2,
        massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                      BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

myModel.ExplicitDynamicsStep(
    name='folding',
    previous=('Initial' if SKIP_PINCHING else 'pinching'),
    timePeriod=3.0,
    linearBulkViscosity=0.06,
    quadBulkViscosity=1.2,
    massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                  BELOW_MIN, 1, 0, 0.0, 0.0, 1, None), ))

try:
    myModel.fieldOutputRequests['F-Output-1'].setValues(
        timeInterval=0.01,
        variables=('CF', 'RF', 'RM', 'RT', 'SE', 'SF', 'U', 'UR', 'UT'))
except KeyError:
    myModel.FieldOutputRequest('F-Output-1',
        createStepName=FIRST_STEP, timeInterval=0.01,
        variables=('CF', 'RF', 'RM', 'RT', 'SE', 'SF', 'U', 'UR', 'UT'))

try:
    myModel.historyOutputRequests['H-Output-1'].setValues(
        timeInterval=0.01,
        variables=('ALLAE', 'ALLCD', 'ALLCW', 'ALLDC', 'ALLDMD', 'ALLFD',
                   'ALLIE', 'ALLKE', 'ALLMW', 'ALLPD', 'ALLPW', 'ALLSE',
                   'ALLVD', 'ALLWK'))
except KeyError:
    myModel.HistoryOutputRequest('H-Output-1',
        createStepName=FIRST_STEP, timeInterval=0.01,
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
# [10] LOAD & BOUNDARY CONDITIONS
# =============================================================
tol_bc = 0.5

# --- Initial BCs on shell nodes ---
v_bc1 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc,       yMin=-R_y - tol_bc, zMin=L_total - tol_bc,
    xMax= tol_bc,       yMax= R_y + tol_bc, zMax=L_total + tol_bc)
a.Set(name='BC1-nodes', vertices=v_bc1)
print('BC-1 (Y-extreme at Z=L): %d vertices' % len(v_bc1))

myModel.DisplacementBC(
    name='BC-1', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC1-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc2 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc,       zMin=L_total - tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc,       zMax=L_total + tol_bc)
a.Set(name='BC2-nodes', vertices=v_bc2)
print('BC-2 (X-extreme at Z=L): %d vertices' % len(v_bc2))

myModel.DisplacementBC(
    name='BC-2', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC2-nodes'],
    u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc3 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc,       yMin=-R_y - tol_bc, zMin=-tol_bc,
    xMax= tol_bc,       yMax= R_y + tol_bc, zMax= tol_bc)
a.Set(name='BC3-nodes', vertices=v_bc3)
print('BC-3 (Y-extreme at Z=0): %d vertices' % len(v_bc3))

myModel.DisplacementBC(
    name='BC-3', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC3-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc4 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc,       zMin=-tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc,       zMax= tol_bc)
a.Set(name='BC4-nodes', vertices=v_bc4)
print('BC-4 (X-extreme at Z=0): %d vertices' % len(v_bc4))

myModel.DisplacementBC(
    name='BC-4', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC4-nodes'],
    u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# --- Viscous pressure (pinching step, propagates to folding) ---
allFaces = inst.faces[:]
myModel.Pressure(
    amplitude=UNSET,
    createStepName=FIRST_STEP,
    distributionType=VISCOUS,
    field='',
    magnitude=1.4623e-06,
    name='Load-1',
    region=Region(side1Faces=allFaces))

# --- Concentrated pinching forces ---
y_threshold_top = 10.0

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

if not SKIP_PINCHING:
    myModel.ConcentratedForce(
        name='pinchload-top',
        createStepName='pinching',
        region=a.sets['pinch-top'],
        cf1=0.0, cf2=-PINCH_FORCE, cf3=0.0,
        amplitude='Pinching')

    myModel.ConcentratedForce(
        name='pinchload-bot',
        createStepName='pinching',
        region=a.sets['pinch-bot'],
        cf1=0.0, cf2=PINCH_FORCE, cf3=0.0,
        amplitude='Pinching')

# --- Folding rotation BCs ---
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

print('===== Loads & BCs applied =====')

# =============================================================
# [11] MESH MODULE
# =============================================================
# Mesh seed sizes are injectable (validation mesh-convergence study) but
# default to the production values, so normal MBO/production runs are
# byte-for-byte unchanged.
try:
    MESH_GLOBAL_SIZE
except NameError:
    MESH_GLOBAL_SIZE = 2      # global part seed (mm)
try:
    MESH_LOCAL_SIZE
except NameError:
    MESH_LOCAL_SIZE = 1.0     # refined hinge-region edge seed (mm)

myPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=MESH_GLOBAL_SIZE)

tol_mesh = 0.5
rect_edges = myPart.edges.getByBoundingBox(
    xMin=-R_x - tol_mesh, yMin=y_min_rect - tol_mesh, zMin=z_min_rect - tol_mesh,
    xMax= R_x + tol_mesh, yMax=y_max_rect + tol_mesh, zMax=z_max_rect + tol_mesh)

if len(rect_edges) > 0:
    myPart.seedEdgeBySize(edges=rect_edges, size=MESH_LOCAL_SIZE,
                          deviationFactor=0.1, constraint=FINER)
    print('Seeded %d edges in bounding rectangle with size=%s'
          % (len(rect_edges), MESH_LOCAL_SIZE))

inner_faces = myPart.faces.getByBoundingBox(
    xMin=-R_x - tol_mesh, yMin=y_min_rect - tol_mesh, zMin=z_min_rect - tol_mesh,
    xMax= R_x + tol_mesh, yMax=y_max_rect + tol_mesh, zMax=z_max_rect + tol_mesh)
inner_face_ids = set(f.index for f in inner_faces)
outer_faces = [f for f in myPart.faces if f.index not in inner_face_ids]

if len(inner_faces) > 0:
    myPart.setMeshControls(regions=inner_faces, elemShape=QUAD_DOMINATED,
                           technique=FREE)
    print('Inner faces: %d -> QUAD_DOMINATED FREE' % len(inner_faces))

structured_cnt = 0
sweep_cnt = 0
for face in outer_faces:
    region = (face, )
    try:
        myPart.setMeshControls(regions=region, elemShape=QUAD,
                               technique=STRUCTURED)
        structured_cnt += 1
    except:
        myPart.setMeshControls(regions=region, elemShape=QUAD,
                               technique=SWEEP,
                               algorithm=ADVANCING_FRONT)
        sweep_cnt += 1

print('Outer faces: %d STRUCTURED, %d SWEEP' % (structured_cnt, sweep_cnt))

myPart.setElementType(
    elemTypes=(
        ElemType(elemCode=S3R, elemLibrary=EXPLICIT,
                 secondOrderAccuracy=ON, hourglassControl=STIFFNESS),
    ),
    regions=(inner_faces, ))

myPart.setElementType(
    elemTypes=(
        ElemType(elemCode=S4R, elemLibrary=EXPLICIT,
                 secondOrderAccuracy=ON, hourglassControl=STIFFNESS),
    ),
    regions=(tuple(outer_faces), ))

myPart.generateMesh()

print('Nodes: %d' % len(myPart.nodes))
print('Elements: %d' % len(myPart.elements))
print('===== Mesh generated =====')

a.regenerate()

# =============================================================
# [12] JOB MODULE — output to pinching/sim_%d/
# =============================================================
if not os.path.exists(SIM_DIR):
    os.makedirs(SIM_DIR)
    print('Created directory: %s' % SIM_DIR)

os.chdir(SIM_DIR)
print('Working directory: %s' % os.getcwd())

job_name = 'cutout_points_%d' % csv_num
mdb.Job(
    activateLoadBalancing=False,
    atTime=None,
    contactPrint=OFF,
    description='TapeSpring Hinge - Pinching/Folding',
    echoPrint=OFF,
    explicitPrecision=DOUBLE_PLUS_PACK,
    historyPrint=OFF,
    memory=90,
    memoryUnits=PERCENTAGE,
    model='Model-1',
    modelPrint=OFF,
    multiprocessingMode=DEFAULT,
    name=job_name,
    nodalOutputPrecision=FULL,
    numCpus=4,
    numDomains=4,
    parallelizationMethodExplicit=DOMAIN,
    queue=None,
    resultsFormat=ODB,
    scratch='',
    type=ANALYSIS,
    userSubroutine='',
    waitHours=0,
    waitMinutes=0)

mdb.jobs[job_name].writeInput()
print('===== Job [%s] INP written =====' % job_name)
print('===== All files -> %s =====' % SIM_DIR)
print('===== Script completed successfully! =====')
