#! /user/bin/python
#-*-coding: UTF-8-*-
# -*- coding: mbcs -*-
# =============================================================
#  MODULE 2: PROPERTY
#  - General Stiffness Section, Section Assignment, Material Orientation
#  - Prerequisite: 01_part.py
#  - Run: abaqus cae noGUI=02_property.py
# =============================================================
from abaqus import *
from abaqusConstants import *
from section import *
from material import *
import regionToolset

myModel = mdb.models['Model-1']
myPart  = myModel.parts['TapeSpring_3D']

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
