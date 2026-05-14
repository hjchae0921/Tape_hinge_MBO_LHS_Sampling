#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  MODULE 5: MESH + JOB
#  - Mesh seeding, element types, mesh generation, job creation
#  - Prerequisite: 01_part.py ~ 04_step_interaction_load.py
#  - Run: abaqus cae noGUI=05_mesh.py
# =============================================================
from abaqus import *
from abaqusConstants import *
from mesh import *
from job import *
import regionToolset
import os
import math

L_total = 220.0
D       = 38.0
R       = D / 2.0
csv_num = 1

try:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE_DIR = os.getcwd()

myModel = mdb.models['Model-1']
myPart  = myModel.parts['TapeSpring_3D']
a       = myModel.rootAssembly

# Conservative upper bounds for bounding box
R_x = 25.0
R_y = 25.0

# Recompute bounding rect from cutout data
# (re-read from part edges or use stored values)
# For standalone module, recalculate from part bounding box of cutout region.
# We use part-level edge queries instead.

# =============================================================
# [11] MESH MODULE
# =============================================================
# Global seed
myPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=2)

# Fine seeding in bounding rectangle region
# Re-derive bounding rect from the cutout points CSV
import csv as csv_mod
csv_file = os.path.join(BASE_DIR, 'specimen', 'cutout_points_%d.csv' % csv_num)
points_2d = []
with open(csv_file, 'r') as f:
    reader = csv_mod.reader(f)
    for i, row in enumerate(reader):
        if i == 0 and row[0].strip().upper() == 'X':
            continue
        if len(row) < 2:
            continue
        aa, bb = row[0].strip(), row[1].strip()
        if aa and bb:
            points_2d.append((float(aa), float(bb)))

N_pts = len(points_2d) // 4
top_right_block = points_2d[0:N_pts]
top_left_block  = points_2d[N_pts:2*N_pts]
bot_left_block  = points_2d[2*N_pts:3*N_pts]
bot_right_block = points_2d[3*N_pts:4*N_pts]
top_2d = top_left_block + top_right_block
bot_2d = bot_left_block[::-1] + bot_right_block

all_pts_uv = []
for x2d, y2d in (top_2d + bot_2d):
    u = y2d
    v = x2d + L_total / 2.0
    all_pts_uv.append((u, v))

z_min_cut = min(v for (u, v) in all_pts_uv)
z_max_cut = max(v for (u, v) in all_pts_uv)
u_max_cut = max(abs(u) for (u, v) in all_pts_uv)

rect_spacing = 1.0
z_min_rect = z_min_cut - rect_spacing
z_max_rect = z_max_cut + rect_spacing
y_min_rect = -(u_max_cut + rect_spacing)
y_max_rect = u_max_cut + rect_spacing

tol = 0.5
rect_edges = myPart.edges.getByBoundingBox(
    xMin=-R_x - tol, yMin=y_min_rect - tol, zMin=z_min_rect - tol,
    xMax= R_x + tol, yMax=y_max_rect + tol, zMax=z_max_rect + tol)

if len(rect_edges) > 0:
    myPart.seedEdgeBySize(edges=rect_edges, size=1.0,
                          deviationFactor=0.1, constraint=FINER)
    print('Seeded %d edges in bounding rectangle with size=1.0' % len(rect_edges))

# Face classification: inner (bounding rect) vs outer
inner_faces = myPart.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=y_min_rect - tol, zMin=z_min_rect - tol,
    xMax= R_x + tol, yMax=y_max_rect + tol, zMax=z_max_rect + tol)
inner_face_ids = set(f.index for f in inner_faces)
outer_faces = [f for f in myPart.faces if f.index not in inner_face_ids]

if len(inner_faces) > 0:
    myPart.setMeshControls(regions=inner_faces, elemShape=QUAD_DOMINATED,
                           technique=FREE)
    print('Inner faces (bounding rect): %d -> QUAD_DOMINATED FREE' % len(inner_faces))

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

# Element type: inner S3R, outer S4R
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
# [12] JOB MODULE
# =============================================================
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

print('===== Job created (Double Precision + Full Nodal Output) =====')
print('===== Script completed successfully! =====')
