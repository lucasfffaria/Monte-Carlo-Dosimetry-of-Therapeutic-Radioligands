# -*- coding: utf-8 -*-
"""Article 1, Sections 3.1-3.5 (Tables 1, 4 and 5; Figs. 3b,c, 5 and 6). Writes article1_results.json.
Article 1 analyses: spongiosa-method effective dose (ICRP 116 convention), Ka, H*(10)/Ka,
time-dependent H*(10) for 177Lu-PSMA vs measurements, time-integrated doses for the four nuclides."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postproc_core import *
OUTJ=os.path.join(os.path.dirname(os.path.abspath(__file__)),'article1_results.json')
LAM={n:np.log(2)/D.NUCLIDE_DATA[n]['thalf_h'] for n in NUCS}
R={}
# ------------------------------------------------------------------ targets, spongiosa method
def target_organs_s(run,key):
    out=target_organs(run,key)
    T=tallies(run); lab=META[run]['tally_index']['F906']
    sex='AM' if 'AM' in key else 'AF'
    mk='AM_M' if sex=='AM' else 'AM_F'; tk='TM50_M' if sex=='AM' else 'TM50_F'
    V,E=arr(T[906]); org=MASS[key]['organs']; acc={}
    for i,l in enumerate(lab):
        tg,s=l.split(':'); s=int(s); m=org[str(s)]['mass_g']; f=MeV/(m*1e-3)*1e12
        d=V[i]*f; e=V[i]*E[i]*f
        for name,kk in (('RBMs',mk),('BSEs',tk)):
            w=D.SKEL_MASS[s][kk] if s in D.SKEL_MASS else 0.0
            a=acc.setdefault((tg,name),[0,0,0.0]); a[0]=a[0]+w*d; a[1]=a[1]+(w*e)**2; a[2]+=w
    for (tg,name),(a,s2,ws) in acc.items(): out[tg][name]=(a/ws,np.sqrt(s2)/ws)
    return out
def E_parts(am,af,pos,method):
    dm=dict(am[pos+'AM']); df=dict(af[pos+'AF'])
    if method=='s':
        for d in (dm,df): d['RBM']=d['RBMs']; d['BSE']=d['BSEs']
    return effective_parts(dm,df)
def resp(V,S,nuc):
    a=AN[nuc]; return V@a, np.sqrt((S**2)@(a**2))
w4=weights_from_model(lu_model(4,0.0)); w4n=w4/w4.sum()
RESP={}     # (quantity,position) -> nuc -> (R_r, err_r) per decay, pSv or pGy
TRUNS=((30,'T01_AM_g030','T05_AF_g030'),(50,'T02_AM_g050','T06_AF_g050'),(100,'T03_AM_g100','T07_AF_g100'),(200,'T04_AM_g200','T08_AF_g200'))
print('Effective dose rate, uSv/h per GBq in body (4 h distribution): spongiosa method [DRF method]')
print('  gap  pos   '+''.join(f'{n:>22s}' for n in NUCS))
for gap,ra,rf in TRUNS:
    am=target_organs_s(ra,'TAMstd'); af=target_organs_s(rf,'TAFstd')
    for pos in ('FRONT','SIDE','BACK'):
        row=''
        for meth in ('s','d'):
            Ev,Es=E_parts(am,af,pos,meth)
            for n in NUCS: RESP[(f'E{meth}',pos,gap,n)]=resp(Ev,Es,n)
        for n in NUCS:
            Rs,Ss=RESP[('Es',pos,gap,n)]; Rd,_=RESP[('Ed',pos,gap,n)]
            v=w4n@Rs*GBQH; s=np.sqrt((w4n**2)@(Ss**2))*GBQH; vd=w4n@Rd*GBQH
            R.setdefault('E',{})[f'{pos}_{gap}_{n}']=dict(spong=v,relerr=s/v,drf=vd)
            row+=f'  {v:8.3f}({s/v*100:3.1f}%)[{vd/v:5.3f}]'
        print(f'  {gap:3d}  {pos:5s}'+row)
# ------------------------------------------------------------------ P0: H*(10) and Ka
P=tallies('P0_paciente_em_pe'); labH=META['P0_paciente_em_pe']['tally_index']['F304']; labK=labH
VH,EH=arr(P[304]); VK,EK=arr(P[314]); assert labH==labK
for i,l in enumerate(labH):
    for n in NUCS:
        RESP[('H',l,n)]=resp(VH[i],VH[i]*EH[i],n); RESP[('K',l,n)]=resp(VK[i],VK[i]*EK[i],n)
print('\nKa (uGy/h per GBq) and H*(10)/Ka (Sv/Gy), 4 h distribution')
for l in ('front_g10_z110','front_g30_z110','front_g100_z110','front_g300_z110','left_g100_z110','back_g100_z110'):
    row=''
    for n in NUCS:
        h=w4n@RESP[('H',l,n)][0]*GBQH; k=w4n@RESP[('K',l,n)][0]*GBQH
        R.setdefault('Ka',{})[f'{l}_{n}']=dict(Ka=k,H=h,ratio=h/k); row+=f'   {k:8.3f} H/K {h/k:5.3f}'
    print(f'  {l:18s}'+row)
# bare point-source coefficients at 1 m (no attenuation)
def coef(e,X,Y): return float(np.exp(np.interp(np.log(e),np.log(X),np.log(Y))))
GAM={}
for n in NUCS:
    h=sum(y*coef(e,D.HSTAR_E,D.HSTAR_F) for e,y in D.nuclide_lines(n))/(4*np.pi*1e4)*GBQH
    k=sum(y*coef(e,D.KA_E,D.KA_F) for e,y in D.nuclide_lines(n))/(4*np.pi*1e4)*GBQH
    GAM[n]=(h,k); R.setdefault('Gamma',{})[n]=dict(Hstar=h,Ka=k)
print('\nBare point source at 1 m: Gamma_H*(10) [uSv m2/(h GBq)], Gamma_Ka [uGy m2/(h GBq)], photons/decay')
for n in NUCS: print(f'  {n}: {GAM[n][0]:7.3f}  {GAM[n][1]:7.3f}   yield {sum(y for e,y in D.nuclide_lines(n)):.4f}')
# ------------------------------------------------------------------ time model (177Lu-PSMA-617 typical; other nuclides: same biology)
def wb_dc(t): return lu_model(t)['WB']*np.exp(LAM['Lu177']*t)
def ecum(t): return 0.0 if t<=0 else max(0.0,100.0-wb_dc(t))
def weights_t(t,tv,nuc,bladder=True):
    B=np.exp(-LAM['Lu177']*t)*(ecum(t)-ecum(tv)) if bladder else 0.0
    m=lu_model(max(t,0.5),B)
    if m['REM']<0: m['REM']=0.0
    return weights_from_model(m)*np.exp(-(LAM[nuc]-LAM['Lu177'])*t)
def rate(key,nuc,t,tv,bladder=True):
    return weights_t(t,tv,nuc,bladder)@RESP[key][0]*GBQH      # per GBq administered
# H*(10) at 1 m and 2 m vs time, 7.4 GBq 177Lu; bladder empty vs full (2-h voiding from 0)
print('\nH*(10) (uSv/h), 7.4 GBq 177Lu-PSMA, front z=110: 1 m [empty | before void, 2-h voids] ; 2 m [empty | full] ; 0.5 m')
demir1={0:48,1:38,2:33,4:23,6:15,18:11,24:7,48:5,120:1}; demir2={0:24,1:21,2:17,4:11,6:8,18:5,24:3,48:1,120:0.3}
TC=[]
for t in (0.5,1,2,4,6,18,24,48,72,120,168):
    tv=2*np.floor((t-1e-9)/2) if t>2 else 0.0
    h1e=rate(('H','front_g100_z110','Lu177'),'Lu177',t,t)*7.4; h1f=rate(('H','front_g100_z110','Lu177'),'Lu177',t,tv)*7.4
    h2e=rate(('H','front_g200_z110','Lu177'),'Lu177',t,t)*7.4; h2f=rate(('H','front_g200_z110','Lu177'),'Lu177',t,tv)*7.4
    h05=rate(('H','front_g50_z110','Lu177'),'Lu177',t,t)*7.4
    TC.append(dict(t=t,H1_empty=h1e,H1_full=h1f,H2_empty=h2e,H2_full=h2f,H05_empty=h05))
    print(f'  t={t:6.1f} h  1 m {h1e:6.2f} | {h1f:6.2f}   2 m {h2e:6.2f} | {h2f:6.2f}   0.5 m {h05:6.2f}   Demir 1 m {demir1.get(int(t),"-")} 2 m {demir2.get(int(t),"-")}')
R['timecurve']=TC
# ------------------------------------------------------------------ time-integrated doses from release
A0={'Lu177':7.4,'Tb161':7.4,'Pb212':0.2,'Ac225':0.008}
def integral(key,nuc,t_rel,dv=4.0,bladder=True,t_end=720.0):
    tot=0.0; t0=t_rel
    while t0<t_end:
        t1=min(t0+dv,t_end); ts=np.linspace(t0,t1,41 if t0<96 else 9)
        ys=[rate(key,nuc,t,t0,bladder) for t in ts]; tot+=np.trapezoid(ys,ts); t0=t1
    return tot          # uSv (or uGy) per GBq administered
print('\nTime-integrated doses from release (t_rel) to 30 d, continuous presence (occupancy 1); 4-h voiding; per typical activity')
INT={}
for n in NUCS:
    for t_rel in (4.0,24.0):
        Hf=integral(('H','front_g100_z110',n),n,t_rel)*A0[n]
        Ef=integral(('Es','FRONT',100,n),n,t_rel)*A0[n]; Es=integral(('Es','SIDE',100,n),n,t_rel)*A0[n]
        Eb=integral(('Es','BACK',100,n),n,t_rel)*A0[n]; E50=integral(('Es','FRONT',50,n),n,t_rel)*A0[n]
        E30=integral(('Es','FRONT',30,n),n,t_rel)*A0[n]
        Efe=integral(('Es','FRONT',100,n),n,t_rel,bladder=False)*A0[n]
        # point source, physical decay only (NRC-type), activity at release, r = 1 m from the source
        Hpt=GAM[n][0]*A0[n]*np.exp(-LAM[n]*t_rel)/LAM[n]
        INT[f'{n}_{t_rel:g}']=dict(A0=A0[n],H1m=Hf,E1m_front=Ef,E1m_side=Es,E1m_back=Eb,E50_front=E50,E30_front=E30,E1m_front_nobladder=Efe,H_point_phys=Hpt)
        print(f'  {n} {A0[n]:6.3f} GBq, t_rel {t_rel:4.0f} h: H*(10) 1 m {Hf:9.1f} uSv | E 1 m front {Ef:8.1f} side {Es:7.1f} back {Eb:8.1f}'
              f' | E front 0.5 m {E50:8.1f} 0.3 m {E30:8.1f} | no-bladder E1m {Efe:7.1f} | point/phys H* {Hpt:9.1f}')
R['integrals']=INT
# fraction of the integral delivered in the first day / week after release (177Lu, E front 1 m)
for n in NUCS:
    tot=integral(('Es','FRONT',100,n),n,4.0); d1=integral(('Es','FRONT',100,n),n,4.0,t_end=28.0); w1=integral(('Es','FRONT',100,n),n,4.0,t_end=172.0)
    print(f'  {n}: fraction of E(1 m front, from 4 h) in first 24 h {d1/tot:5.3f}, first 7 d {w1/tot:5.3f}')
    R.setdefault('fractions',{})[n]=dict(day1=d1/tot,week1=w1/tot)
# release criterion ANSN 3.05: H*(10) < 30 uSv/h at 2 m -> max activity in body (4 h distribution), all faces
print('\nActivity in body (GBq) giving 30 uSv/h at 2 m (4 h distribution), front/back/left; and at typical activity')
for n in NUCS:
    hs=[w4n@RESP[('H',f'{f}_g200_z110',n)][0]*GBQH for f in ('front','back','left','right')]
    R.setdefault('release2m',{})[n]=dict(H2m_perGBq=max(hs),Amax=30/max(hs))
    print(f'  {n}: H*(10) 2 m per GBq front/back/left/right ' + ' '.join(f'{h:6.3f}' for h in hs) + f'  -> A_max = {30/max(hs):7.2f} GBq')
json.dump(R,open(OUTJ,'w',newline='\n'),indent=1,default=float)
