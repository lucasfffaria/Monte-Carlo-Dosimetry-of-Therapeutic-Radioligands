# -*- coding: utf-8 -*-
"""Article 1, Sections 2.4 and 3.1-3.3 (Table 2; Fig. 3b,c; Table 3; Fig. 4): H*(10) of the model patient versus
time (7.4 GBq [177Lu]Lu-PSMA-617) for the benchmark, H*(10) per GBq in the body around the patient, comparison with
the simplified sources (in-body point and line sources; bare point source without attenuation), and Tables 2 and 3."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postproc_core import *
run='P0_paciente_em_pe'; T=tallies(run); lab=META[run]['tally_index']['F304']
V,E=arr(T[304]); S=V*E      # pSv per photon [det,region,E]
Vk,Ek=arr(T[314]); Sk=Vk*Ek
def det(name): return lab.index(name)
IA=7.4  # GBq
print('177Lu-PSMA (typical model), H*(10) at 1 m in front of the patient (z=110 cm), 7.4 GBq:')
for t,bl,cmt in ((0.5,0.0,'bladder empty'),(4,0.0,'just after voiding'),(4,11.4,'before void (2-h voids)'),(6,0.0,'after voiding'),(6,9.2,'before void')):
    m=lu_model(t,bl); w=weights_from_model(m)
    val,err=fold(V[det('front_g100_z110')],S[det('front_g100_z110')],'Lu177',w)
    print(f'  t={t:4.1f} h ({cmt:24s}) body {m["WB"]:5.1f} %IA + bladder {bl:4.1f}%: H*(10) = {val*GBQH*IA:6.1f} +- {err*GBQH*IA:4.1f} uSv/h')
print('  Demir et al. 2016 (7.4 GBq, n=23): 48+-13 (start), 23+-6 (4 h), 15+-4 (6 h) uSv/h at 1 m')
# --- source-model comparison per GBq at 4 h distribution (activity normalised to 1 GBq in body)
m=lu_model(4,0.0); w=weights_from_model(m); w=w/w.sum()
print('\nH*(10) per GBq in the body (uSv/h/GBq), 177Lu-PSMA distribution at 4 h, vs simplified point/line sources:')
print('  detector             '+''.join(f'{n:>22s}' for n in NUCS))
for d in ('front_g10_z110','front_g30_z110','front_g50_z110','front_g100_z110','front_g200_z110','front_g300_z110','left_g100_z110','back_g100_z110','front_g100_z160','front_g100_z45'):
    i=det(d); row=''
    for n in NUCS:
        vb,eb=fold(V[i],S[i],n,w); vp,_=fold(V[i],S[i],n,w_point()); vl,_=fold(V[i],S[i],n,w_line())
        row+=f'  {vb*GBQH:8.2f}({eb/vb*100:3.1f}%) P/V {vp/vb:4.2f} L/V {vl/vb:4.2f}'[:22].rjust(22) if False else f'  {vb*GBQH:7.2f} P{vp/vb:4.2f} L{vl/vb:4.2f}'
    print(f'  {d:20s}'+row)
# relative errors at 3 m for realistic distribution
for d in ('front_g300_z110','left_g300_z110','diag45_R225_z110'):
    i=det(d); print('  rel.err',d,' '.join(f'{n} {fold(V[i],S[i],n,w)[1]/fold(V[i],S[i],n,w)[0]*100:.1f}%' for n in NUCS))

# ---------------------------------------------------------------- bare point source (inverse-square law, no attenuation)
print()
def hstar_coef(e):
    return float(np.exp(np.interp(np.log(e),np.log(D.HSTAR_E),np.log(D.HSTAR_F))))
src=np.array([0.0,0.0,105.0])     # PNT source position (global): trunk centre, z=105
P=tallies('P0_paciente_em_pe'); lab=META['P0_paciente_em_pe']['tally_index']['F304']; V,E=arr(P[304])
dets={d[0]:np.array(d[1]) for d in META['P0_paciente_em_pe']['detectors']}
w=weights_from_model(lu_model(4,0.0)); w=w/w.sum()
print('H*(10) per GBq (uSv/h): voxel patient (4 h distribution) vs bare point source in air at the trunk centre')
print('  detector            dist(cm) '+''.join(f'{n:>21s}' for n in NUCS))
for dname in ('front_g10_z110','front_g30_z110','front_g50_z110','front_g100_z110','front_g200_z110','front_g300_z110','left_g100_z110','back_g100_z110'):
    i=lab.index(dname); r=np.linalg.norm(dets[dname]-src); row=''
    for n in NUCS:
        bare=sum(y*hstar_coef(e) for e,y in D.nuclide_lines(n))/(4*np.pi*r*r)*GBQH
        vb,_=fold(V[i],V[i]*E[i],n,w)
        row+=f'   {vb*GBQH:7.2f} bare/vox {bare/(vb*GBQH):4.2f}'
    print(f'  {dname:18s} {r:7.1f}'+row)

# ================================================================ article 1 tables reproduced from the model
# Table 2: time-integrated activity coefficients (TIAC, h) and %IA at 4 h of the reference biokinetic model
print('\nTable 2 (article 1): TIAC (h) and %IA at 4 h of the reference [177Lu]Lu-PSMA-617 model (physical decay included)')
_t = np.concatenate([np.linspace(0, 10, 100001), np.linspace(10, 6000, 600001)[1:]])
def _comp(m):
    return {'Whole body (excl. bladder contents)': m['WB'], 'Kidneys': m['KID'], 'Salivary glands': m['SAL'],
            'Liver': m['LIV'], 'Spleen': m['SPL'], 'Red marrow': m['RM'], 'Tumour (bone, lymph nodes, prostate)': m['BONE'] + m['LN'] + m['PRO'],
            'Remainder (incl. lacrimal glands)': m['REM']}
_mt, _m4 = _comp(lu_model(_t)), _comp(lu_model(4.0))
for k in _mt:
    print(f'  {k:38s} TIAC {np.trapezoid(_mt[k], _t) / 100.0:7.2f} h   %IA(4 h) {_m4[k]:6.2f}')
# Table 3: H*(10) per GBq in the body at 110 cm height, by direction and gap from the body surface
print('\nTable 3 (article 1): H*(10) per GBq in the body (uSv/h/GBq), 4 h distribution, bladder empty, z = 110 cm')
print('  radionuclide   direction   10 cm      30 cm      50 cm        1 m        2 m        3 m')
_V3, _E3 = arr(tallies('P0_paciente_em_pe')[304]); _w3 = weights_from_model(lu_model(4, 0.0)); _w3 = _w3 / _w3.sum()
for n in NUCS:
    for face, nome in (('front', 'Front'), ('back', 'Back'), ('left', 'Left')):
        vals = [fold(_V3[det(f'{face}_g{g}_z110')], (_V3 * _E3)[det(f'{face}_g{g}_z110')], n, _w3)[0] * GBQH for g in (10, 30, 50, 100, 200, 300)]
        print(f'  {n:13s}  {nome:9s}' + ''.join(f'{v:11.3g}' for v in vals))

# ================================================================ derived datasets behind Table 3 and Figs. 4, S2b and S3
_HERE = os.path.dirname(os.path.abspath(__file__))
# (1) every detector of the patient-alone run: H*(10) (H) and air kerma (K) per GBq in the body (4 h distribution),
#     relative error of H, and H*(10) per GBq for the in-body point (Hp) and line (Hl) sources; xyz = sphere centre (cm)
_P = tallies('P0_paciente_em_pe'); _lab = META['P0_paciente_em_pe']['tally_index']['F304']
_VH, _EH = arr(_P[304]); _VK, _EK = arr(_P[314]); _dets = {d[0]: d[1] for d in META['P0_paciente_em_pe']['detectors']}
_out = {}
for i, l in enumerate(_lab):
    row = {'xyz': _dets[l]}
    for n in NUCS:
        a = AN[n]; Rr = _VH[i] @ a; Sr = np.sqrt(((_VH[i] * _EH[i]) ** 2) @ (a ** 2))
        row[n] = dict(H=float(_w3 @ Rr * GBQH), rel=float(np.sqrt((_w3 ** 2) @ (Sr ** 2)) / (_w3 @ Rr)),
                      Hp=float(w_point() @ Rr * GBQH), Hl=float(w_line() @ Rr * GBQH), K=float(_w3 @ (_VK[i] @ a) * GBQH))
    _out[l] = row
json.dump(_out, open(os.path.join(_HERE, 'article1_detector_results.json'), 'w'))
print(f'\nwritten: article1_detector_results.json ({len(_out)} detectors; uSv/h or uGy/h per GBq in the body)')
# (2) Fig. S2b: photon fluence rate per unit lethargy (cm-2 s-1 GBq-1) in the spheres 1 m from the body surface at 110 cm
#     height (tally F324: front, left and back), [177Lu]Lu-PSMA-617 distribution at 4 h
_t = _P[324]; _d = _t['dims']
_v = _t['vals'][:, 0].reshape([_d[k] for k in 'fdusmcet'])[:, 0, :NE, 0, 0, 0, :50, :NREG] / (PE[None, :, None, None] * (1.0 / NREG))
_edges = np.concatenate([[0.010], np.geomspace(0.012, 2.7, 50)]); _leth = np.log(_edges[1:] / _edges[:-1])
_cols = {f'{n}_{face}': np.einsum('kbr,r,k->b', _v[c], _w3, AN[n]) / _leth * 1e9 for n in NUCS for c, face in enumerate(('front', 'left', 'back'))}
with open(os.path.join(_HERE, 'article1_photon_spectra_1m.csv'), 'w', encoding='utf-8', newline='\n') as f:
    f.write('E_low_MeV,E_high_MeV,E_centre_keV,' + ','.join(_cols) + '\n')
    for b in range(50):
        f.write(f'{_edges[b]:.5g},{_edges[b + 1]:.5g},{1000 * np.sqrt(_edges[b] * _edges[b + 1]):.5g},' + ','.join(f'{_cols[k][b]:.5g}' for k in _cols) + '\n')
print('written: article1_photon_spectra_1m.csv (fluence rate per unit lethargy, cm-2 s-1 GBq-1)')
