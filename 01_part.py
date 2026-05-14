#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  MODULE 1: PART
#  - Parameters, CSV read, base shell, cutout, partitions
#  - Run: abaqus cae noGUI=01_part.py
# =============================================================
from abaqus import *
from abaqusConstants import *
import __main__
from part import *
from sketch import *
import regionToolset
import os
import math
import csv

# =============================================================
# [0] Parameters
# =============================================================
L_total = 220.0
D       = 38.0
R       = D / 2.0
csv_num = 10

try:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE_DIR = os.getcwd()
csv_file = os.path.join(BASE_DIR, 'specimen', 'cutout_points_%d.csv' % csv_num)

COUPLING_STRIP = 5.0     # mm from each end for coupling surface
PINCH_OFFSET   = 5.0     # mm offset for pinching grid

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
# [2] PART MODULE - Base Elliptical Shell
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
print('Base elliptical shell created (Rx=%.1f, Ry=%.1f, L=%.1f)' % (R_x, R_y, L_total))

# =============================================================
# [3] 2D -> 3D projection
# =============================================================
N = len(points_2d) // 4
top_right_block = points_2d[0:N]
top_left_block  = points_2d[N:2*N]
bot_left_block  = points_2d[2*N:3*N]
bot_right_block = points_2d[3*N:4*N]

top_2d = top_left_block + top_right_block
bot_2d = bot_left_block[::-1] + bot_right_block

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
# [4] Datum Plane (YZ) + sketch + CutExtrude (both +/-X)
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
# --- Basic partitions (XY/YZ/XZ at origin) ---
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

# --- Cutout bounding rectangle partitions ---
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
print('Bounding rect (spacing=%.1f): Z=[%.1f, %.1f], Y=[%.1f, %.1f]' %
      (rect_spacing, z_min_rect, z_max_rect, y_min_rect, y_max_rect))

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

# --- NEW: Coupling strip partitions (5 mm from each end) ---
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

# --- NEW: Pinching grid partitions ---
#   XY planes at Z = L/2-5, L/2, L/2+5
#   YZ planes at X = +5, -5
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
