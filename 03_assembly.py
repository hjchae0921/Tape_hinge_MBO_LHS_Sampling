#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  MODULE 3: ASSEMBLY
#  - Instance, Reference Points, Surface Coupling, Point Inertia
#  - Prerequisite: 01_part.py, 02_property.py
#  - Run: abaqus cae noGUI=03_assembly.py
# =============================================================
from abaqus import *
from abaqusConstants import *
from assembly import *
import regionToolset
import math

L_total = 220.0
D       = 38.0
R       = D / 2.0
COUPLING_STRIP = 5.0

myModel = mdb.models['Model-1']
myPart  = myModel.parts['TapeSpring_3D']

# Recompute R_x, R_y (needed for face selection tolerances)
# Read eccentricity from the first datum csys or re-derive from part
# For robustness, use bounding box of part
bb = myPart.queryGeometry()
# R_x and R_y can be inferred, but we just use a generous tolerance
R_x = 25.0   # conservative upper bound
R_y = 25.0

# =============================================================
# [7] ASSEMBLY MODULE
# =============================================================
a = myModel.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(dependent=ON, name='TapeSpring_3D-1',
           part=myPart)

print('===== Assembly created =====')

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
# Coupling surface: face strip within COUPLING_STRIP mm from each end
tol = 0.1

# Top coupling surface (Z = L-5 ~ L)
facesA = inst.faces.getByBoundingBox(
    xMin=-R_x - tol, yMin=-R_y - tol, zMin=L_total - COUPLING_STRIP - tol,
    xMax= R_x + tol, yMax= R_y + tol, zMax=L_total + tol)

# Bottom coupling surface (Z = 0 ~ 5)
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

print('===== Surface Coupling created =====')

# --- Point Inertia on RP-C ---
a.engineeringFeatures.PointMassInertia(
    name='Inertia-1',
    region=a.sets['RP-C'],
    mass=1e-06,
    i11=1e-06, i22=1e-06, i33=1e-06)

print('===== Point Inertia on RP-C created =====')
print('===== Assembly module completed =====')
