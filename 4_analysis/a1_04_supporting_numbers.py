# -*- coding: utf-8 -*-
"""Article 1: supporting numbers quoted in the text (height and direction dependence of the field, Section 3.2;
ratios of the in-body sources and of the bare point source to the voxel E, Sections 3.3-3.4; release criterion per
unit H*(10), Section 3.5; relative uncertainties). Requires article1_results.json (run a1_03 first)."""
import os, sys, json, numpy as np
_HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,_HERE)
_SRC=open(os.path.join(_HERE,'a1_03_article1_results.py'),encoding='utf-8').read()
_PRE,_POS=_SRC.split('# ------------------------------------------------------------------ time model')
exec(_PRE.replace("print(","(lambda *a,**k: None)("))                        # response matrices (printing suppressed)
exec('# ---- time model'+_POS.split("# H*(10) at 1 m and 2 m vs time")[0])      # time model (rate, weights_t)

# ================================================================ height/direction profiles and source-model ratios
lab=labH
def h(l,n,w): return w@RESP[('H',l,n)][0]*GBQH
print('Height dependence, front g100 (uSv/h per GBq, 4 h): z=20..160')
for n in NUCS: print(' ',n,' '.join(f'{h(f"front_g100_z{z}",n,w4n):6.2f}' for z in (20,45,70,90,110,135,160)))
print('left g100 heights 20..160:'); 
for n in NUCS: print(' ',n,' '.join(f'{h(f"left_g100_z{z}",n,w4n):6.2f}' for z in (20,45,70,90,110,135,160)))
print('Direction at z=110: front/back/left/right g100 and diag R125 45/135/225/315')
for n in NUCS: print(' ',n,' '.join(f'{h(l,n,w4n):6.2f}' for l in ('front_g100_z110','back_g100_z110','left_g100_z110','right_g100_z110','diag45_R125_z110','diag135_R125_z110','diag225_R125_z110','diag315_R125_z110')))
print('Per-GBq-in-body H*(10) and E vs distribution time (no bladder):')
for t in (0.5,4,24,72):
    m=lu_model(t,0.0); w=weights_from_model(m); wn=w/w.sum()
    row=f'  t={t:5.1f}h  H10 {h("front_g10_z110","Lu177",wn):6.2f}  H30 {h("front_g30_z110","Lu177",wn):6.2f}  H100 {h("front_g100_z110","Lu177",wn):5.2f}  H300 {h("front_g300_z110","Lu177",wn):5.3f}'
    row+=f'  E30F {wn@RESP[("Es","FRONT",30,"Lu177")][0]*GBQH:6.3f}  E100F {wn@RESP[("Es","FRONT",100,"Lu177")][0]*GBQH:6.3f}  E30B {wn@RESP[("Es","BACK",30,"Lu177")][0]*GBQH:6.3f}'
    row+=f'  | Tb H30 {h("front_g30_z110","Tb161",wn):6.2f} Pb H30 {h("front_g30_z110","Pb212",wn):7.1f}'
    print(row)
wp=w_point(); wl=w_line()
print('E ratios in-body point/voxel and line/voxel (spongiosa), by position and gap')
for n in NUCS:
    for pos in ('FRONT','SIDE','BACK'):
        r=[]
        for gap in (30,50,100,200):
            Rs=RESP[('Es',pos,gap,n)][0]; v=w4n@Rs; r.append(f'{(wp@Rs)/v:5.2f}/{(wl@Rs)/v:5.2f}')
        print(f'  {n} {pos:5s} '+'  '.join(r))
# bare point E-estimate: bare point H*(10) at the exposed person's front surface vs voxel E
print('Bare-point H*(10) at the gap point (inverse square from trunk centre) / voxel E')
for n in NUCS:
    for pos,face in (('FRONT','front'),('SIDE','left'),('BACK','back')):
        r=[]
        for gap in (30,50,100,200):
            l=f'{face}_g{gap}_z110'; i=lab.index(l)
            c=np.array([d[1] for d in META['P0_paciente_em_pe']['detectors'] if d[0]==l][0]); dist=np.linalg.norm(c-np.array([0,0,105.]))
            bare=GAM[n][0]*(100.0/dist)**2
            r.append(f'{bare/(w4n@RESP[("Es",pos,gap,n)][0]*GBQH):5.2f}')
        print(f'  {n} {pos:5s} '+' '.join(r))

# ================================================================ release criterion per unit H*(10)
print()
A0={'Lu177':7.4,'Tb161':7.4,'Pb212':0.2,'Ac225':0.008}
J=json.load(open(os.path.join(_HERE,'article1_results.json'),encoding='utf-8'))
for n in NUCS:
    h2=rate(('H','front_g200_z110',n),n,4.0,4.0)*A0[n]; e1=rate(('Es','FRONT',100,n),n,4.0,4.0)*A0[n]
    Einf=J['integrals'][f'{n}_4']['E1m_front']; Eb=J['integrals'][f'{n}_4']['E1m_front_nobladder']
    print(f'{n}: H*(10) 2 m front at 4 h {h2:7.3f} uSv/h ; E rate 1 m front at 4 h {e1:7.3f} uSv/h ; Einf {Einf:7.2f} uSv -> Einf/H2m {Einf/h2:6.1f} h ; Einf/E1m-rate {Einf/e1:6.1f} h ; bladder share {1-Eb/Einf:5.3f}')

# ================================================================ relative uncertainties of H*(10)
print()
for n in NUCS:
    row=[]
    for l in ('front_g10_z110','front_g30_z110','front_g100_z110','back_g100_z110','left_g100_z110','front_g300_z110','left_g300_z110','back_g300_z110'):
        R_,S_=RESP[('H',l,n)]; v=w4n@R_; s=np.sqrt((w4n**2)@(S_**2)); row.append(f'{l.split("_z")[0]}:{s/v*100:4.1f}%')
    print(n,' '.join(row))

# ================================================================ Fig. 5a: H*(10) at the survey-meter point / E of the exposed person
print('\nFig. 5a (article 1): H*(10) at the survey-meter point (patient-alone run) / sex-averaged E of the exposed person')
print('  position        gap     ' + '  '.join(f'{n:>6s}' for n in NUCS))
for pos, face, name in (('FRONT', 'front', 'Face to face'), ('BACK', 'back', 'Behind'), ('SIDE', 'left', 'Side by side')):
    for gap in (30, 50, 100, 200):
        r = [(w4n @ RESP[('H', f'{face}_g{gap}_z110', n)][0]) / (w4n @ RESP[('Es', pos, gap, n)][0]) for n in NUCS]
        print(f'  {name:13s} {gap:4d} cm  ' + '  '.join(f'{x:6.2f}' for x in r))
# ================================================================ Fig. 6b: cumulative share of the 30-day E 1 m in front (release at 4 h)
_cum = {}
for n in NUCS:
    key = ('Es', 'FRONT', 100, n); T = []; Y = []; t0 = 4.0
    while t0 < 720.0:
        t1 = min(t0 + 4.0, 720.0); ts = np.linspace(t0, t1, 41 if t0 < 96 else 9)
        T.append(ts); Y.append(np.array([rate(key, n, t, t0) for t in ts])); t0 = t1
    T = np.concatenate(T); Y = np.concatenate(Y)
    c = np.concatenate([[0], np.cumsum(0.5 * (Y[1:] + Y[:-1]) * np.diff(T))]); _cum[n] = (T - 4.0, c / c[-1])
with open(os.path.join(_HERE, 'article1_cumulative_fraction.csv'), 'w', encoding='utf-8', newline='\n') as f:
    f.write('t_after_release_h,' + ','.join(NUCS) + '\n')
    for k in range(len(_cum['Lu177'][0])):
        f.write(f'{_cum["Lu177"][0][k]:.6g},' + ','.join(f'{_cum[n][1][k]:.6g}' for n in NUCS) + '\n')
print('\nwritten: article1_cumulative_fraction.csv (Fig. 6b; E 1 m in front, release at 4 h, voiding every 4 h)')
