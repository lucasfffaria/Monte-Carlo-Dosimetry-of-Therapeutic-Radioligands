# -*- coding: utf-8 -*-
"""Article 1, Section 3.1 and Fig. 3a: verification of the dosimetric chain against ICRP 116.
Ratio of the computed to the reference E/Phi for broad parallel AP, PA, LLAT and RLAT beams (runs V01-V08), with the
spongiosa approximation (convention of the reference coefficients) and with the skeletal dose-response functions (DRF)."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postproc_core import *
def vdoses(run,sex):
    T=tallies(run); key='TAMstd' if sex=='AM' else 'TAFstd'; area=META[run]['beam_area']; nE=len(META[run]['energies'])
    tags=[t[0] for t in D.TISSUES]
    def a(t):
        d=t['dims']; v=t['vals'][:,0].reshape(d['f'],d['u']); e=t['vals'][:,1].reshape(d['f'],d['u'])
        return v[:,:nE],e[:,:nE]
    v,e=a(T[6]); dose={}; err={}
    for i,tag in enumerate(tags):
        mg=MASS[key]['tissues'][tag]['mass_g']; dose[tag]=v[i]*MeV/(mg*1e-3)*1e12*area*nE; err[tag]=e[i]
    mk='AM_M' if sex=='AM' else 'AM_F'; tk='TM50_M' if sex=='AM' else 'TM50_F'
    rb=[s for s in D.SPONGIOSA if (s,'AM') in D.DRF]; tm=[s for s in sorted(set(D.SPONGIOSA+D.MEDULLARY)) if (s,'TM50') in D.DRF]
    for name,sites,base,kk in (('RBM',rb,104,mk),('BSE',tm,404,tk)):
        s_=0; ws=0
        for q,s in enumerate(sites):
            vv,_=a(T[base+10*q]); w=D.SKEL_MASS[s][kk]; s_=s_+w*vv[0]*area*nE; ws+=w
        dose[name]=s_/ws
    # spongiosa method (ICRP 116 reference for RBM/endosteum): mass-weighted spongiosa absorbed dose
    vs,_=a(T[906]); org=MASS[key]['organs']
    sp=[s for s in D.SPONGIOSA if str(s) in org]
    for name,kk in (('RBMs',mk),('BSEs',tk)):
        acc=0; ws=0
        for j,s in enumerate(sp):
            m=org[str(s)]['mass_g']; w=D.SKEL_MASS[s][kk] if s in D.SKEL_MASS else 0
            acc=acc+w*vs[j]*MeV/(m*1e-3)*1e12*area*nE; ws+=w
        dose[name]=acc/ws
    return dose,err
def Eff(dm,df,rbm='RBM',bse='BSE'):
    E=0
    for t,w in WT.items():
        if t=='REM': h=np.mean([0.5*(dm[x]+df[x]) for x in REM],axis=0)
        elif t=='RBM': h=0.5*(dm[rbm]+df[rbm])
        elif t=='BSE': h=0.5*(dm[bse]+df[bse])
        else: h=0.5*(dm[t]+df[t])
        E=E+w*h
    return E
ref=D.EPHI     # ICRP 116, Photons_Effective dose.txt (E/Phi, pSv cm2): columns E, AP, PA, LLAT, RLAT, ROT, ISO
cols={'AP':1,'PA':2,'LLAT':3,'RLAT':4}
En=META['V01_AM_AP']['energies']
pairs={'AP':('V01_AM_AP','V05_AF_AP'),'PA':('V02_AM_PA','V06_AF_PA'),'LLAT':('V03_AM_LLAT','V07_AF_LLAT'),'RLAT':('V04_AM_RLAT','V08_AF_RLAT')}
R={}
for g,(ra,rf) in pairs.items():
    dm,em=vdoses(ra,'AM'); df,ef=vdoses(rf,'AF')
    E1=Eff(dm,df); E2=Eff(dm,df,'RBMs','BSEs'); r=np.interp(En,ref[:,0],ref[:,cols[g]])
    R[g]=(E1/r,E2/r)
print('E/Phi ratio MCNP/ICRP116 (DRF method | spongiosa method), production V runs')
print('  E(MeV) '+'  '.join(f'{g:>13s}' for g in pairs))
for i,en in enumerate(En):
    print(f'  {en:6.3f} '+'  '.join(f'{R[g][0][i]:6.3f}|{R[g][1][i]:5.3f}' for g in pairs))
allr=np.array([R[g][1] for g in pairs]); sel=np.array(En)>=0.03
print('spongiosa method, 30 keV-3 MeV: mean %.3f, min %.3f, max %.3f'%(allr[:,sel].mean(),allr[:,sel].min(),allr[:,sel].max()))
allr=np.array([R[g][0] for g in pairs])
print('DRF method,       30 keV-3 MeV: mean %.3f, min %.3f, max %.3f'%(allr[:,sel].mean(),allr[:,sel].min(),allr[:,sel].max()))
