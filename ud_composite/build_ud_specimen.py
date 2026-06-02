#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  TapeSpring Hinge UD — Build (Pinching + Folding)
#
#  Geometry / BC / load / mesh: identical to tube_hinge_pinching.py.
#  Material: CU250 UD (LAMINA + FailStress) + CompositeShellSection.
#  Output: ud_composite/sim_L<layup_id>_S<csv_num>/
#
#  Runner injects:
#    csv_num    : int (specimen, e.g. 197 or 291)
#    layup_id   : int (1..74, from ud_layups.LAYUPS)
#    BASE_DIR   : repo root (parent of ud_composite/)
#
#  Run via runner — `abaqus cae noGUI=_temp_build_L<id>_S<num>.py`.
# =============================================================
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
import regionToolset
import sys, os, math, csv

# =============================================================
# [0] Parameters
# =============================================================
L_total = 220.0
D       = 38.0
R       = D / 2.0

try:
    csv_num
except NameError:
    csv_num = 291
try:
    layup_id
except NameError:
    layup_id = 11   # default sanity: [45/-45]

COUPLING_STRIP = 5.0
PINCH_OFFSET   = 5.0
PINCH_FORCE    = 3.0
FOLD_ANGLE     = 1.4835

# Resolve BASE_DIR (repo root, parent of ud_composite/)
try:
    BASE_DIR
except NameError:
    try:
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        BASE_DIR = os.getcwd()

UD_DIR       = os.path.join(BASE_DIR, 'ud_composite')
SPECIMEN_DIR = os.path.join(BASE_DIR, 'specimen')

# Import layup table (ud_composite/ must be on sys.path)
if UD_DIR not in sys.path:
    sys.path.insert(0, UD_DIR)
from ud_layups import get_layup, PLY_THICKNESS

LAY = get_layup(layup_id)
PLY_ANGLES = LAY['angles']
N_PLIES    = len(PLY_ANGLES)
TOT_THICK  = N_PLIES * PLY_THICKNESS

SIM_DIR  = os.path.join(UD_DIR, 'sim_L%d_S%d' % (layup_id, csv_num))
csv_file = os.path.join(SPECIMEN_DIR, 'cutout_points_%d.csv' % csv_num)

print('============================================')
print('  UD BUILD  layup_id=%d (%s)  specimen=%d' % (layup_id, LAY['label'], csv_num))
print('  plies=%d  total_t=%.3f mm' % (N_PLIES, TOT_THICK))
print('  angles=%s' % str(PLY_ANGLES))
print('  SIM_DIR=%s' % SIM_DIR)
print('============================================')

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
# [3] 2D -> 3D projection (identical to tube_hinge_pinching.py)
# =============================================================
if len(points_2d) >= 2 and points_2d[0] == points_2d[-1]:
    points_2d = points_2d[:-1]

zero_idx = [i for i in range(len(points_2d)) if abs(points_2d[i][1]) < 1e-6]
if len(zero_idx) >= 2:
    i_left  = zero_idx[0]
    i_right = zero_idx[1]
    top_2d = points_2d[i_left : i_right + 1]
    bot_2d_ccw = points_2d[i_right:] + [points_2d[i_left]]
    bot_2d = list(reversed(bot_2d_ccw))
else:
    top_2d = [p for p in points_2d if p[1] >= 0]
    bot_2d = [p for p in points_2d if p[1] <= 0]

top_pts = [(y2d, x2d + L_total / 2.0) for (x2d, y2d) in top_2d]
bot_pts = [(y2d, x2d + L_total / 2.0) for (x2d, y2d) in bot_2d]
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
# [4] CutExtrude (both +X and -X)
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
        return
    left_curve  = sk_pts[:best_start + 1]
    right_curve = sk_pts[best_end:]
    if len(left_curve) >= 2:
        sketch.Spline(points=tuple(left_curve))
    sketch.Line(point1=sk_pts[best_start], point2=sk_pts[best_end])
    if len(right_curve) >= 2:
        sketch.Spline(points=tuple(right_curve))


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

s_cut = myModel.ConstrainedSketch(name='__cutout__', sheetSize=500.0, transform=t)
myPart.projectReferencesOntoSketch(sketch=s_cut, filter=COPLANAR_EDGES)

top_sk = [(v - L_total / 2.0, u) for (u, v) in top_pts]
bot_sk = [(v - L_total / 2.0, u) for (u, v) in bot_pts]
if len(top_sk) > 0 and len(bot_sk) > 0:
    bot_sk[0]  = top_sk[-1]
    bot_sk[-1] = top_sk[0]

draw_half_curve(s_cut, top_sk)
draw_half_curve(s_cut, bot_sk)

myPart.CutExtrude(
    sketchPlane=myPart.datums[dp_id],
    sketchUpEdge=myPart.datums[da_id],
    sketchPlaneSide=SIDE1,
    sketchOrientation=TOP,
    sketch=s_cut, flipExtrudeDirection=ON,
    depth=R_x + 5.0)
del myModel.sketches['__cutout__']

s_cut2 = myModel.ConstrainedSketch(name='__cutout2__', sheetSize=500.0, transform=t)
myPart.projectReferencesOntoSketch(sketch=s_cut2, filter=COPLANAR_EDGES)
draw_half_curve(s_cut2, top_sk)
draw_half_curve(s_cut2, bot_sk)
myPart.CutExtrude(
    sketchPlane=myPart.datums[dp_id],
    sketchUpEdge=myPart.datums[da_id],
    sketchPlaneSide=SIDE1,
    sketchOrientation=TOP,
    sketch=s_cut2, flipExtrudeDirection=OFF,
    depth=R_x + 5.0)
del myModel.sketches['__cutout2__']
print('Cutout extrudes applied (+X and -X).')

# =============================================================
# [5] Partitions
# =============================================================
dp_xy = myPart.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
dp_yz = myPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
dp_xz = myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)

for dp_obj, name in [(dp_xy, 'XY'), (dp_yz, 'YZ'), (dp_xz, 'XZ')]:
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp_obj.id],
            faces=myPart.faces[:])
    except Exception:
        pass

all_cut_pts = top_pts + bot_pts
z_min_cut = min(v for (u, v) in all_cut_pts)
z_max_cut = max(v for (u, v) in all_cut_pts)
u_max_cut = max(abs(u) for (u, v) in all_cut_pts)

rect_spacing = 1.0
z_min_rect = z_min_cut - rect_spacing
z_max_rect = z_max_cut + rect_spacing
y_min_rect = -(u_max_cut + rect_spacing)
y_max_rect =  (u_max_cut + rect_spacing)

for offset_val, plane_type in [
    (z_min_rect, XYPLANE), (z_max_rect, XYPLANE),
    (y_min_rect, XZPLANE), (y_max_rect, XZPLANE),
]:
    dp = myPart.DatumPlaneByPrincipalPlane(principalPlane=plane_type, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id], faces=myPart.faces[:])
    except Exception:
        pass

for offset_val in (COUPLING_STRIP, L_total - COUPLING_STRIP):
    dp = myPart.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id], faces=myPart.faces[:])
    except Exception:
        pass

z_mid = L_total / 2.0
for offset_val, plane_type in [
    (z_mid - PINCH_OFFSET, XYPLANE),
    (z_mid,                XYPLANE),
    (z_mid + PINCH_OFFSET, XYPLANE),
    ( PINCH_OFFSET,        YZPLANE),
    (-PINCH_OFFSET,        YZPLANE),
]:
    dp = myPart.DatumPlaneByPrincipalPlane(principalPlane=plane_type, offset=offset_val)
    try:
        myPart.PartitionFaceByDatumPlane(
            datumPlane=myPart.datums[dp.id], faces=myPart.faces[:])
    except Exception:
        pass

print('Part [TapeSpring_3D] created.')

# =============================================================
# [6] PROPERTY MODULE — CU250 UD + CompositeShellSection
# =============================================================
# Consistent units: tonne, mm, N, s, MPa.
#   E [MPa = N/mm2], rho [t/mm3] = kg/m3 * 1e-12
#   1984 kg/m3 = 1.984e-9 t/mm3
RHO_UD = 1.984e-9     # t/mm3
E1     = 135000.0     # MPa (135 GPa)
E2     = 8600.0       # MPa
NU12   = 0.30
G12    = 4100.0       # MPa
G13    = 2000.0       # MPa
G23    = 2000.0       # MPa

# Strength [MPa]
XT = 2916.0
XC = 1127.0
YT = 55.0
YC = 180.0
S_SHEAR = 113.0
# FailStress table: Xt, Xc, Yt, Yc, S, cross-product limit Fb*, biaxial-stress limit
# Abaqus FailStress signature:
#   (sigma_T1, sigma_C1, sigma_T2, sigma_C2, sigma_12, biaxial f*, biaxial sigma_lim)
# Default 0.0 disables the cross-product limit (=> standard Tsai-Wu F12 = -0.5*sqrt(F11*F22)).
FAIL_STRESS = ((XT, XC, YT, YC, S_SHEAR, 0.0, 0.0),)

mat = myModel.Material(name='CU250_UD')
mat.Density(table=((RHO_UD,),))
mat.Elastic(type=LAMINA, table=((E1, E2, NU12, G12, G13, G23),))
mat.elastic.FailStress(table=FAIL_STRESS)

# Composite Shell Section -- one layer per ply
section_layers = []
for ply_idx, theta in enumerate(PLY_ANGLES):
    section_layers.append(
        SectionLayer(material='CU250_UD',
                     thickness=PLY_THICKNESS,
                     orientAngle=float(theta),
                     numIntPts=3,
                     plyName='Ply_%d' % (ply_idx + 1))
    )

myModel.CompositeShellSection(
    name='UD_Section',
    preIntegrate=OFF,
    idealization=NO_IDEALIZATION,
    thicknessType=UNIFORM,
    poissonDefinition=DEFAULT,
    temperature=GRADIENT,
    useDensity=OFF,
    integrationRule=SIMPSON,
    layup=section_layers)

myPart.SectionAssignment(
    offset=0.0,
    offsetField='',
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=myPart.faces[:]),
    sectionName='UD_Section',
    thicknessAssignment=FROM_SECTION)

# Material orientation: cylindrical CSYS; ply 1-direction = global Z (axial).
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

print('Property: CU250_UD LAMINA + %d-ply CompositeShellSection.' % N_PLIES)

# =============================================================
# [7] ASSEMBLY MODULE (identical to original)
# =============================================================
a = myModel.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(dependent=ON, name='TapeSpring_3D-1', part=myModel.parts['TapeSpring_3D'])
inst = a.instances['TapeSpring_3D-1']

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

tol = 0.1
facesA = inst.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=-R_y - tol, zMin=L_total - COUPLING_STRIP - tol,
    xMax= R_x + tol, yMax= R_y + tol, zMax=L_total + tol)
facesB = inst.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=-R_y - tol, zMin=-tol,
    xMax= R_x + tol, yMax= R_y + tol, zMax=COUPLING_STRIP + tol)
a.Surface(name='CoupSurf-Top', side1Faces=facesA)
a.Surface(name='CoupSurf-Bot', side1Faces=facesB)

myModel.Coupling(
    controlPoint=Region(referencePoints=(rpA,)),
    couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None,
    name='Constraint-1', surface=a.surfaces['CoupSurf-Top'],
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
myModel.Coupling(
    controlPoint=Region(referencePoints=(rpB,)),
    couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None,
    name='Constraint-2', surface=a.surfaces['CoupSurf-Bot'],
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a.engineeringFeatures.PointMassInertia(
    name='Inertia-1', region=a.sets['RP-C'],
    mass=1e-06, i11=1e-06, i22=1e-06, i33=1e-06)
print('Assembly + couplings done.')

# =============================================================
# [8] STEP MODULE (pinching 1 s + folding 3 s)
# =============================================================
myModel.SmoothStepAmplitude(
    name='Pinching', timeSpan=TOTAL,
    data=((0.0, 0.0), (1.0, 1.0), (4.0, 1.0), (5.0, 0.0), (8.0, 0.0)))
myModel.SmoothStepAmplitude(
    name='Rotation', timeSpan=TOTAL,
    data=((0.0, 0.0), (1.0, 0.0), (4.0, 1.0), (5.0, 1.0), (8.0, 0.0)))

myModel.ExplicitDynamicsStep(
    name='pinching', previous='Initial',
    timePeriod=1.0,
    linearBulkViscosity=0.06, quadBulkViscosity=1.2,
    massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                  BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))
myModel.ExplicitDynamicsStep(
    name='folding', previous='pinching',
    timePeriod=3.0,
    linearBulkViscosity=0.06, quadBulkViscosity=1.2,
    massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                  BELOW_MIN, 1, 0, 0.0, 0.0, 1, None), ))

# Field output: include S + CFAILURE for composite failure post-processing.
# (CFAILURE is computed at each section point when FailStress is defined.)
F_VARS = ('S', 'SF', 'SM', 'SE', 'CFAILURE',
          'U', 'UR', 'CF', 'RF', 'RM', 'RT')
try:
    myModel.fieldOutputRequests['F-Output-1'].setValues(
        timeInterval=0.01, variables=F_VARS)
except KeyError:
    myModel.FieldOutputRequest('F-Output-1',
        createStepName='pinching', timeInterval=0.01, variables=F_VARS)

H_VARS = ('ALLAE', 'ALLCD', 'ALLCW', 'ALLDC', 'ALLDMD', 'ALLFD',
          'ALLIE', 'ALLKE', 'ALLMW', 'ALLPD', 'ALLPW', 'ALLSE',
          'ALLVD', 'ALLWK')
try:
    myModel.historyOutputRequests['H-Output-1'].setValues(
        timeInterval=0.01, variables=H_VARS)
except KeyError:
    myModel.HistoryOutputRequest('H-Output-1',
        createStepName='pinching', timeInterval=0.01, variables=H_VARS)

print('Steps + output requests created.')

# =============================================================
# [9] INTERACTION MODULE
# =============================================================
myModel.ContactProperty('IntProp-1')
myModel.interactionProperties['IntProp-1'].TangentialBehavior(formulation=FRICTIONLESS)
myModel.interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, pressureOverclosure=HARD)

myModel.ContactExp(createStepName='Initial', name='Int-1')
myModel.interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
myModel.interactions['Int-1'].contactPropertyAssignments.appendInStep(
    assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

# =============================================================
# [10] LOAD & BOUNDARY CONDITIONS
# =============================================================
tol_bc = 0.5

v_bc1 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc,       yMin=-R_y - tol_bc, zMin=L_total - tol_bc,
    xMax= tol_bc,       yMax= R_y + tol_bc, zMax=L_total + tol_bc)
a.Set(name='BC1-nodes', vertices=v_bc1)
myModel.DisplacementBC(
    name='BC-1', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC1-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc2 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc,       zMin=L_total - tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc,       zMax=L_total + tol_bc)
a.Set(name='BC2-nodes', vertices=v_bc2)
myModel.DisplacementBC(
    name='BC-2', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC2-nodes'],
    u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc3 = inst.vertices.getByBoundingBox(
    xMin=-tol_bc,       yMin=-R_y - tol_bc, zMin=-tol_bc,
    xMax= tol_bc,       yMax= R_y + tol_bc, zMax= tol_bc)
a.Set(name='BC3-nodes', vertices=v_bc3)
myModel.DisplacementBC(
    name='BC-3', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC3-nodes'],
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

v_bc4 = inst.vertices.getByBoundingBox(
    xMin=-R_x - tol_bc, yMin=-tol_bc,       zMin=-tol_bc,
    xMax= R_x + tol_bc, yMax= tol_bc,       zMax= tol_bc)
a.Set(name='BC4-nodes', vertices=v_bc4)
myModel.DisplacementBC(
    name='BC-4', createStepName='Initial',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['BC4-nodes'],
    u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

allFaces = inst.faces[:]
myModel.Pressure(
    amplitude=UNSET,
    createStepName='pinching',
    distributionType=VISCOUS,
    field='',
    magnitude=1.4623e-06,
    name='Load-1',
    region=Region(side1Faces=allFaces))

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

myModel.ConcentratedForce(
    name='pinchload-top', createStepName='pinching',
    region=a.sets['pinch-top'],
    cf1=0.0, cf2=-PINCH_FORCE, cf3=0.0, amplitude='Pinching')
myModel.ConcentratedForce(
    name='pinchload-bot', createStepName='pinching',
    region=a.sets['pinch-bot'],
    cf1=0.0, cf2=PINCH_FORCE, cf3=0.0, amplitude='Pinching')

myModel.DisplacementBC(
    name='BC-5', createStepName='folding',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['RP-A'],
    u1=UNSET, u2=UNSET, u3=UNSET,
    ur1=FOLD_ANGLE, ur2=UNSET, ur3=UNSET, amplitude='Rotation')
myModel.DisplacementBC(
    name='BC-6', createStepName='folding',
    distributionType=UNIFORM, fieldName='', localCsys=None,
    region=a.sets['RP-B'],
    u1=UNSET, u2=UNSET, u3=UNSET,
    ur1=-FOLD_ANGLE, ur2=UNSET, ur3=UNSET, amplitude='Rotation')

# =============================================================
# [11] MESH MODULE
# =============================================================
myPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=2)
tol_mesh = 0.5

rect_edges = myPart.edges.getByBoundingBox(
    xMin=-R_x - tol_mesh, yMin=y_min_rect - tol_mesh, zMin=z_min_rect - tol_mesh,
    xMax= R_x + tol_mesh, yMax=y_max_rect + tol_mesh, zMax=z_max_rect + tol_mesh)
if len(rect_edges) > 0:
    myPart.seedEdgeBySize(edges=rect_edges, size=1.0,
                          deviationFactor=0.1, constraint=FINER)

inner_faces = myPart.faces.getByBoundingBox(
    xMin=-R_x - tol_mesh, yMin=y_min_rect - tol_mesh, zMin=z_min_rect - tol_mesh,
    xMax= R_x + tol_mesh, yMax=y_max_rect + tol_mesh, zMax=z_max_rect + tol_mesh)
inner_face_ids = set(f.index for f in inner_faces)
outer_faces = [f for f in myPart.faces if f.index not in inner_face_ids]

if len(inner_faces) > 0:
    myPart.setMeshControls(regions=inner_faces, elemShape=QUAD_DOMINATED, technique=FREE)

for face in outer_faces:
    region = (face, )
    try:
        myPart.setMeshControls(regions=region, elemShape=QUAD, technique=STRUCTURED)
    except Exception:
        myPart.setMeshControls(regions=region, elemShape=QUAD,
                               technique=SWEEP, algorithm=ADVANCING_FRONT)

myPart.setElementType(
    elemTypes=(ElemType(elemCode=S3R, elemLibrary=EXPLICIT,
                        secondOrderAccuracy=ON, hourglassControl=STIFFNESS), ),
    regions=(inner_faces, ))
myPart.setElementType(
    elemTypes=(ElemType(elemCode=S4R, elemLibrary=EXPLICIT,
                        secondOrderAccuracy=ON, hourglassControl=STIFFNESS), ),
    regions=(tuple(outer_faces), ))

myPart.generateMesh()
print('Mesh: nodes=%d  elems=%d' % (len(myPart.nodes), len(myPart.elements)))
a.regenerate()

# =============================================================
# [12] JOB MODULE — write INP to ud_composite/sim_L<id>_S<num>/
# =============================================================
if not os.path.exists(SIM_DIR):
    os.makedirs(SIM_DIR)
os.chdir(SIM_DIR)

job_name = 'ud_L%d_S%d' % (layup_id, csv_num)
mdb.Job(
    activateLoadBalancing=False, atTime=None,
    contactPrint=OFF, description='UD CompositeShell TapeSpring',
    echoPrint=OFF,
    explicitPrecision=DOUBLE_PLUS_PACK,
    historyPrint=OFF,
    memory=90, memoryUnits=PERCENTAGE,
    model='Model-1', modelPrint=OFF,
    multiprocessingMode=DEFAULT,
    name=job_name,
    nodalOutputPrecision=FULL,
    numCpus=4, numDomains=4,
    parallelizationMethodExplicit=DOMAIN,
    queue=None, resultsFormat=ODB, scratch='',
    type=ANALYSIS, userSubroutine='',
    waitHours=0, waitMinutes=0)

mdb.jobs[job_name].writeInput()
print('===== Job [%s] INP written -> %s =====' % (job_name, SIM_DIR))
print('===== UD build complete =====')
