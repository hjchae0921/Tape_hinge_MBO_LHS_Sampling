# -*- coding: utf-8 -*-
# Seed the new (folded-FI) campaign's initial LHS rows 0..11 from the EXISTING
# production ODBs (root sim_<n>), WITHOUT re-solving. Computes MAX_FI_FOLDED
# (last folding frame), MAX_FI_ALL/PINCHING/FOLDING, MAX_SE; writes to
# new_simulation/results.csv.  Also serves as a smoke-test of the folded post.
#
# Injected globals: ROOT_DIR, NEW_DIR, SPECS (list)
from __future__ import print_function
import os, csv
from abaqus import session
from abaqusConstants import *
from math import sqrt as _sqrt

try: ROOT_DIR
except NameError: ROOT_DIR = os.getcwd()
try: NEW_DIR
except NameError: NEW_DIR = os.path.join(ROOT_DIR, 'new_simulation')
try: SPECS
except NameError: SPECS = list(range(0, 12))

SPECDIR = os.path.join(NEW_DIR, 'specimen')
RESULTS = os.path.join(NEW_DIR, 'results.csv')

sF1t=139.47; sF1c=63.42; sF3=17.73; sF4=5.07
sJ1=1.0/sF1t-1.0/sF1c; sK11=1.0/(sF1t*sF1c); sK33=1.0/(sF3**2); sK12=-sK11/2.0
HEADER=['SPECIMEN_NUM','MAX_FI_ALL','MAX_FI_PINCHING','MAX_FI_FOLDING',
        'MAX_FI_FOLDED','MAX_SE','Eccentricity','y1','x1','Dx','y3','yP2','fP3',
        'x2_derived','x3_derived','xP3_derived','yP3_derived']


def _root(N, rp, rm):
    c=[r for r in (rp,rm) if (r>0 if N>=0 else r<0)]
    if not c: return None
    if len(c)==1: return c[0]
    return c[0] if abs(c[0])>=abs(c[1]) else c[1]


def frame_max_fi(frame, inst):
    av=list(frame.fieldOutputs.keys())
    if 'SF' not in av or 'SM' not in av: return None
    fV=frame.fieldOutputs['SF'].getSubset(region=inst).values
    mV=frame.fieldOutputs['SM'].getSubset(region=inst).values
    mMap={}
    for mv in mV: mMap[(int(mv.elementLabel),int(mv.integrationPoint))]=mv.data
    fmax=-1e99
    for fv in fV:
        k=(int(fv.elementLabel),int(fv.integrationPoint))
        if k not in mMap: continue
        N1,N2,N3,N12,N13,N23=fv.data; M1,M2,M12=mMap[k]
        Nx=0.5*(N1+N2)+N12; Ny=0.5*(N1+N2)-N12; Nxy=0.5*(N2-N1)
        Mx=0.5*(M1+M2)+M12; My=0.5*(M1+M2)-M12; M=Mx if abs(Mx)>=abs(My) else My
        fIP=sJ1*(Nx+Ny)+sK11*(Nx**2+Ny**2)+sK12*Nx*Ny+sK33*Nxy**2
        if fIP>=1.0:
            fC=1.0+abs(M)/sF4
        else:
            dx=(sJ1+sK12*Ny)**2-4*sK11*(sJ1*Ny+sK11*Ny**2+sK33*Nxy**2-1)
            sFx=_root(Nx,(-(sJ1+sK12*Ny)+_sqrt(dx))/(2*sK11),(-(sJ1+sK12*Ny)-_sqrt(dx))/(2*sK11)) if dx>=0 else None
            fCx=(fIP+abs(M)/sF4) if (sFx is None or abs(sFx)<1e-9) else abs(Nx/sFx)+abs(M)/sF4
            dy=(sJ1+sK12*Nx)**2-4*sK11*(sJ1*Nx+sK11*Nx**2+sK33*Nxy**2-1)
            sFy=_root(Ny,(-(sJ1+sK12*Nx)+_sqrt(dy))/(2*sK11),(-(sJ1+sK12*Nx)-_sqrt(dy))/(2*sK11)) if dy>=0 else None
            fCy=(fIP+abs(M)/sF4) if (sFy is None or abs(sFy)<1e-9) else abs(Ny/sFy)+abs(M)/sF4
            fC=fCx if fCx>=fCy else fCy
        if fC>fmax: fmax=fC
    return fmax if fmax>-1e99 else None


def read_dv(spec):
    p=os.path.join(SPECDIR,'cutout_points_%d.csv'%spec)
    d={'Eccentricity':'','y1':'','x1':'','Dx':'','y3':'','yP2':'','fP3':'',
       'x2_derived':'','x3_derived':'','xP3_derived':'','yP3_derived':''}
    rows=list(csv.reader(open(p,'r')))
    if len(rows)>1:
        dr=rows[1]; keys=['Eccentricity','y1','x1','Dx','y3','yP2','fP3','x2_derived','x3_derived','xP3_derived','yP3_derived']
        for i,k in enumerate(keys, start=2):
            if len(dr)>i and dr[i].strip(): d[k]=float(dr[i].strip())
    return d


out_rows=[]
for spec in SPECS:
    odbp=os.path.join(ROOT_DIR,'sim_%d'%spec,'cutout_points_%d.odb'%spec)
    if not os.path.isfile(odbp):
        print('missing', spec); continue
    if odbp in session.odbs.keys(): session.odbs[odbp].close()
    o=session.openOdb(name=odbp, readOnly=True)
    inst=o.rootAssembly.instances[[k for k in o.rootAssembly.instances.keys() if k.upper()!='ASSEMBLY'][0]]
    # objective only needs folded FI (last folding frame) + max ALLSE -> fast.
    folded=''; se=[]
    if 'folding' in o.steps:
        st=o.steps['folding']
        folded=frame_max_fi(st.frames[-1], inst)
        if folded is None: folded=''
    for sn in [s for s in ['pinching','folding'] if s in o.steps]:
        st=o.steps[sn]
        for rk in st.historyRegions.keys():
            reg=st.historyRegions[rk]
            if 'ALLSE' in reg.historyOutputs:
                for tv,dv in reg.historyOutputs['ALLSE'].data: se.append(dv)
    o.close()
    dvs=read_dv(spec)
    row={'SPECIMEN_NUM':spec,'MAX_FI_ALL':'','MAX_FI_PINCHING':'',
         'MAX_FI_FOLDING':'','MAX_FI_FOLDED':folded,'MAX_SE':max(se) if se else ''}
    row.update(dvs)
    out_rows.append(row)
    print('SPEC %d  FOLDED=%s  SE=%s'%(spec,folded,max(se) if se else ''))

# merge into results.csv
existing=[]
if os.path.exists(RESULTS):
    existing=list(csv.DictReader(open(RESULTS,'r')))
by={str(r['SPECIMEN_NUM']):r for r in existing}
for r in out_rows: by[str(r['SPECIMEN_NUM'])]=r
merged=sorted(by.values(), key=lambda r:int(r['SPECIMEN_NUM']))
try: f=open(RESULTS,'w',newline='')
except TypeError: f=open(RESULTS,'wb')
try:
    w=csv.DictWriter(f,fieldnames=HEADER,lineterminator='\n',extrasaction='ignore')
    w.writeheader(); w.writerows(merged)
finally: f.close()
print('seeded %d rows -> %s'%(len(out_rows), RESULTS))
