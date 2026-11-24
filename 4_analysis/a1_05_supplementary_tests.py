# -*- coding: utf-8 -*-
"""Article 1, Sections 2.7, 3.1 and 4 (Supplementary Sections S3.3 and S3.4): the two supplementary runs R1 and R2 (4 November 2026).

R1 - folding of the 161Tb emission lines on the energy grid. In the same run, half of the histories sampled the 21
     photon lines of 161Tb (D) and the other half the same lines folded on the 27-energy grid with the weights of the
     post-processing (G), with the [177Lu]Lu-PSMA-617 activity distribution at 4 h. D/G measures the folding error
     independently of the statistics of the production run; G and D are also compared with the folding of the
     production run P0 (P).
R2 - concrete floor. T = tally with a 30 cm concrete floor under the patient; F = part of it due to photons that have
     left the floor (CF cards); U = T - F. T/U is the increase due to the floor (correlated estimate within the run);
     T/P is the same increase relative to the production run P0 without floor (independent runs); U/P should be ~1.

Reads 2_MCNP_results/R1_P0_Tb161_direto.m, R2_P0_piso_concreto.m and P0_paciente_em_pe.m, and
1_MCNP_inputs/metadados_testes.json (source energies and regions of R1). Writes article1_supplementary_tests.csv.
Units: uSv/h (H*(10)) and uGy/h (air kerma) per GBq in the body, distribution at 4 h just after voiding (Table 3).
Uncertainties: one statistical standard deviation.
"""
import json
import os
import re
import sys

import numpy as np

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import postproc_core as C          # noqa: E402
from mctal import read_mctal       # noqa: E402

R1, R2, P0 = 'R1_P0_Tb161_direto', 'R2_P0_piso_concreto', 'P0_paciente_em_pe'
META = json.load(open(C.paths.INPUTS / 'metadados_testes.json', encoding='utf-8'))
LAB = C.META[P0]['tally_index']['F304']
W4 = C.weights_from_model(C.lu_model(4.0, 0.0))
W4 = W4 / W4.sum()
KEY = ['front_g10_z110', 'front_g30_z110', 'front_g50_z110', 'front_g100_z110', 'front_g200_z110', 'front_g300_z110',
       'back_g10_z110', 'back_g100_z110', 'back_g300_z110', 'left_g10_z110', 'left_g100_z110', 'left_g300_z110',
       'diag45_R225_z110', 'diag225_R225_z110']
HEIGHTS = [20, 45, 70, 90, 110, 135, 160]
CSV = []


def cube(t):
    d = t['dims']
    shp = [d[k] for k in 'fdusmcet']
    return t['vals'][:, 0].reshape(shp), t['vals'][:, 1].reshape(shp)


def ratio(a, sa, b, sb):
    """a/b and its uncertainty (a and b independent)"""
    r = a / b
    return r, abs(r) * np.sqrt((sa / a) ** 2 + (sb / b) ** 2)


def production(tn, nuc):
    """folding of the production run P0: value and uncertainty per GBq at every detector"""
    V, E = C.arr(C.tallies(P0)[tn])
    v, s = C.fold(V, V * E, nuc, W4)
    return v * C.GBQH, s * C.GBQH


def where(label):
    """(direction, gap or distance from the trunk axis, height) of a detector label"""
    m = re.match(r'(front|back|left|right)_g(\d+)_z(\d+)', label)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    m = re.match(r'diag\d+_R(\d+)_z(\d+)', label)
    return 'diag', int(m.group(1)), int(m.group(2))


def near_1m(label):
    d, g, _ = where(label)
    return g <= 100 if d != 'diag' else g <= 125


# ================================================================== R1: 161Tb lines versus their folding on the grid
def analyse_R1():
    T = read_mctal(C.paths.mctal(R1))
    print(f'[{R1}] {T["nps"]:,} histories')
    T = T['tallies']
    m = META[R1]
    ent = m['energy_entries']
    p = np.array([e['p_deck'] for e in ent])
    p = p / p.sum()                                    # sampling probability of each entry (MCNP normalises SP6)
    y = np.array([e['photons_per_decay'] for e in ent])
    line = np.array([e['kind'] == 'line' for e in ent])
    fac = y / p                                        # per source photon -> per decay
    nE, nR = len(ent), len(m['regions'])

    def split(tn, ebin=0, tbin=None):
        """D and G per decay [f] with uncertainties, for time bin tbin (None = sum over the regions)"""
        v, e = cube(T[tn])
        assert v.shape[2] == nE + 1 and v.shape[7] == nR + 1, (tn, v.shape)
        it = nR if tbin is None else tbin
        v = v[:, 0, :nE, 0, 0, 0, ebin, it] * fac[None, :]
        s = v * e[:, 0, :nE, 0, 0, 0, ebin, it]
        return (v[:, line].sum(1), np.sqrt((s[:, line] ** 2).sum(1)),
                v[:, ~line].sum(1), np.sqrt((s[:, ~line] ** 2).sum(1)))

    print('\n' + '=' * 110)
    print('R1 - 161Tb: real lines (D) versus the same lines folded on the grid (G), same run; P = folding of P0')
    print('=' * 110)
    out = {}
    for tn, q in ((304, 'H*(10)'), (314, 'Ka')):
        D, sD, G, sG = (x * C.GBQH for x in split(tn))
        P, sP = production(tn, 'Tb161')
        rDG, sDG = ratio(D, sD, G, sG)
        rGP, sGP = ratio(G, sG, P, sP)
        rDP, sDP = ratio(D, sD, P, sP)
        unit = 'uSv/h/GBq' if tn == 304 else 'uGy/h/GBq'
        print(f'\n{q} ({unit})')
        print(f'  {"detector":20s} {"D (lines)":>16s} {"G (grid)":>16s} {"P (production)":>16s} '
              f'{"D/G":>14s} {"G/P":>14s} {"D/P":>14s}')
        for d in KEY:
            i = LAB.index(d)
            print(f'  {d:20s} {D[i]:9.4g} +-{100 * sD[i] / D[i]:4.2f}% {G[i]:9.4g} +-{100 * sG[i] / G[i]:4.2f}% '
                  f'{P[i]:9.4g} +-{100 * sP[i] / P[i]:4.2f}% {rDG[i]:7.4f}+-{sDG[i]:.4f} {rGP[i]:7.4f}+-{sGP[i]:.4f} '
                  f'{rDP[i]:7.4f}+-{sDP[i]:.4f}')
        w = 1 / sDG ** 2
        nr = np.array([near_1m(d) for d in LAB])
        far = np.abs(rDG - 1) > 2 * sDG
        print(f'  76 detectors: D/G {rDG.min():.4f}-{rDG.max():.4f}, inverse-variance mean {(rDG * w).sum() / w.sum():.4f},'
              f' uncertainty {100 * sDG.min():.2f}-{100 * sDG.max():.2f}%; {nr.sum()} detectors within about 1 m '
              f'(gap <= 100 cm, diagonals <= 125 cm): {rDG[nr].min():.4f}-{rDG[nr].max():.4f}')
        print(f'  |D/G - 1| > 2 SD at {far.sum()} detectors ({", ".join(np.array(LAB)[far])}); expected by chance '
              f'{76 * 0.0455:.1f}.  G/P {rGP.min():.3f}-{rGP.max():.3f} (mean {np.mean(rGP):.3f})')
        out[q] = (rDG, sDG, nr, far, rGP)
        for i, d in enumerate(LAB):
            CSV.append(f'R1,{q},Tb161,{d},{D[i]:.6g},{sD[i]:.3g},{G[i]:.6g},{sG[i]:.3g},{P[i]:.6g},{sP[i]:.3g},'
                       f'{rDG[i]:.5f},{sDG[i]:.5f},{rGP[i]:.5f},{sGP[i]:.5f}')
    D, sD, G, sG = split(344)
    r, sr = ratio(D, sD, G, sG)
    print(f'\nAll photons leaving the patient (F344, world air cell): D/G = {r[0]:.5f} +- {sr[0]:.5f}')
    print('  by source region (D/G):')
    row = []
    for k, tag in enumerate(r_[0] for r_ in m['regions']):
        Dk, sDk, Gk, sGk = split(344, tbin=k)
        if Dk[0] > 0 and Gk[0] > 0:
            rk, srk = ratio(Dk[0], sDk[0], Gk[0], sGk[0])
            row.append(f'{tag} {rk:.4f}+-{srk:.4f}')
    for j in range(0, len(row), 5):
        print('    ' + '   '.join(row[j:j + 5]))
    v, e = cube(T[324])
    edges = np.concatenate([[0.0], np.geomspace(0.012, 2.7, 50)])
    print('\nFluence 1 m from the body surface (F324), D/G by energy band of the photon at the detector:')
    for c, d in enumerate(m['tally_index']['F324']):
        vv = v[c, 0, :nE, 0, 0, 0, :50, nR] * fac[:, None]
        ss = vv * e[c, 0, :nE, 0, 0, 0, :50, nR]
        txt = []
        for lo, hi in ((0, 0.040), (0.040, 0.060), (0.060, 0.100), (0.100, 3.0)):
            b = (edges[1:] > lo) & (edges[1:] <= hi)
            rb, srb = ratio(vv[line][:, b].sum(), np.sqrt((ss[line][:, b] ** 2).sum()),
                            vv[~line][:, b].sum(), np.sqrt((ss[~line][:, b] ** 2).sum()))
            txt.append(f'{1000 * lo:.0f}-{1000 * hi:.0f} keV {rb:.3f}+-{srb:.3f}')
        print(f'  {d:16s} ' + '   '.join(txt))
    print('  (the folding moves the birth energy of the Dy K X-rays, 45-54 keV, to the grid energies, so the bands that'
          '\n   contain them differ as expected; the integrated quantities are the ones that matter)')
    rDG, sDG, nr, far, rGP = out['H*(10)']
    print('\nNumbers quoted in Section 3.1: D/G of H*(10) '
          f'{rDG[nr].min():.3f}-{rDG[nr].max():.3f} at the {nr.sum()} detectors within about 1 m and '
          f'{rDG.min():.2f}-{rDG.max():.2f} at all 76, uncertainties {100 * sDG.min():.1f}-{100 * sDG.max():.1f}%, '
          f'{far.sum()} ratios beyond 2 SD; all photons {r[0]:.4f} +- {sr[0]:.4f}; mean G/P {np.mean(rGP):.3f}')


# ================================================================== R2: concrete floor
def analyse_R2():
    T = read_mctal(C.paths.mctal(R2))
    print(f'\n[{R2}] {T["nps"]:,} histories')
    T = T['tallies']
    m = META[R2]
    regs = [r[0] for r in m['regions']]
    assert regs == C.REG
    rp = np.array([m['region_prob'][r] for r in regs])
    pe = np.array(m['energy_prob'])

    def arr_d(t, di):
        """as postproc_core.arr, for D bin di (0 = total, 1 = photons that have left the floor)"""
        v, e = cube(t)
        assert v.shape[1] == 2, 'R2 without CF cards?'
        v, e = v[:, di, :C.NE, 0, 0, 0, 0, :C.NREG], e[:, di, :C.NE, 0, 0, 0, 0, :C.NREG]
        v = v / (pe[None, :, None] * rp[None, None, :])
        return np.transpose(v, (0, 2, 1)), np.transpose(e, (0, 2, 1))

    print('\n' + '=' * 110)
    print('R2 - 30 cm concrete floor: T = with floor; F = photons that have left the floor; U = T - F; P = P0 (no floor)')
    print('=' * 110)
    res = {}
    for tn, q in ((304, 'H*(10)'), (314, 'Ka')):
        VT, ET = arr_d(T[tn], 0)
        VF, EF = arr_d(T[tn], 1)
        for nuc in C.NUCS:
            t, st = C.fold(VT, VT * ET, nuc, W4)
            f, sf = C.fold(VF, VF * EF, nuc, W4)
            u = t - f
            su = np.sqrt(np.clip(st ** 2 - sf ** 2, 0, None))      # U and F treated as independent
            P, sP = production(tn, nuc)
            P, sP = P / C.GBQH, sP / C.GBQH
            frac = f / t
            sfrac = np.sqrt(u ** 2 * sf ** 2 + f ** 2 * su ** 2) / t ** 2
            tu, stu = 1 / (1 - frac), sfrac / (1 - frac) ** 2
            tp, stp = ratio(t, st, P, sP)
            up, sup = ratio(u, su, P, sP)
            res[(q, nuc)] = dict(t=t * C.GBQH, st=st * C.GBQH, P=P * C.GBQH, sP=sP * C.GBQH, tu=tu, stu=stu,
                                 tp=tp, stp=stp, up=up, sup=sup)
            for i, d in enumerate(LAB):
                CSV.append(f'R2,{q},{nuc},{d},{t[i] * C.GBQH:.6g},{st[i] * C.GBQH:.3g},{f[i] * C.GBQH:.6g},'
                           f'{sf[i] * C.GBQH:.3g},{P[i] * C.GBQH:.6g},{sP[i] * C.GBQH:.3g},{tu[i]:.5f},{stu[i]:.5f},'
                           f'{tp[i]:.5f},{stp[i]:.5f},{up[i]:.5f},{sup[i]:.5f}')
    for q in ('H*(10)', 'Ka'):
        print(f'\n{q}: increase due to the floor, T/U (correlated) and T/P (versus the production run); U/P should be ~1')
        print(f'  {"detector":20s}' + ''.join(f'{n + " T/U":>16s}{"T/P":>14s}' for n in C.NUCS))
        for d in KEY + [f'front_g100_z{h}' for h in HEIGHTS if h != 110] + [f'left_g100_z{h}' for h in HEIGHTS if h != 110]:
            i = LAB.index(d)
            print(f'  {d:20s}' + ''.join(f'{res[(q, n)]["tu"][i]:8.4f}+-{res[(q, n)]["stu"][i]:.4f}'
                                         f'{res[(q, n)]["tp"][i]:8.3f}+-{res[(q, n)]["stp"][i]:.3f}' for n in C.NUCS))
        for n in C.NUCS:
            r = res[(q, n)]
            print(f'  {n}: U/P {r["up"].min():.3f}-{r["up"].max():.3f} (mean {np.mean(r["up"]):.3f})')
    r = res[('H*(10)', 'Lu177')]
    print('\nRatio H*(10)(2 m)/H*(10)(1 m) in front of the patient, 177Lu (Supplementary Section S3.4):')
    for h in (90, 110, 135):
        i1, i2 = LAB.index(f'front_g100_z{h}'), LAB.index(f'front_g200_z{h}')
        a, sa = ratio(r['P'][i2], r['sP'][i2], r['P'][i1], r['sP'][i1])
        b, sb = ratio(r['t'][i2], r['st'][i2], r['t'][i1], r['st'][i1])
        print(f'  z = {h:3d} cm: without floor {a:.3f} +- {sa:.3f}   with floor {b:.3f} +- {sb:.3f}')
    print('\nH*(10) 1 m in front of the patient, 177Lu (uSv/h/GBq), by detector height:')
    for h in HEIGHTS:
        i = LAB.index(f'front_g100_z{h}')
        print(f'  z = {h:3d} cm: without floor {r["P"][i]:6.3f} +- {r["sP"][i]:.3f}   with floor {r["t"][i]:6.3f} +- '
              f'{r["st"][i]:.3f}   increase {100 * (r["tu"][i] - 1):5.1f} +- {100 * r["stu"][i]:.1f} %')
    # ranges quoted in Sections 3.1 and 4 (H*(10), 110 cm height unless stated, four radionuclides)
    inc = {n: 100 * (res[('H*(10)', n)]['tu'] - 1) for n in C.NUCS}

    def rng(sel):
        idx = [LAB.index(d) for d in LAB if sel(d)]
        vals = np.concatenate([inc[n][idx] for n in C.NUCS])
        return vals.min(), vals.max()
    z110 = lambda d: where(d)[2] == 110 and where(d)[0] != 'diag'    # noqa: E731
    groups = [('within 50 cm', lambda d: z110(d) and where(d)[1] <= 50),
              ('1 m, front or back', lambda d: z110(d) and where(d)[1] == 100 and where(d)[0] in ('front', 'back')),
              ('1 m, beside (left or right)', lambda d: z110(d) and where(d)[1] == 100 and where(d)[0] in ('left', 'right')),
              ('2 m, front or back', lambda d: z110(d) and where(d)[1] == 200 and where(d)[0] in ('front', 'back')),
              ('2 m, beside (left or right)', lambda d: z110(d) and where(d)[1] == 200 and where(d)[0] in ('left', 'right')),
              ('1 m, all directions', lambda d: z110(d) and where(d)[1] == 100),
              ('2 m, all directions', lambda d: z110(d) and where(d)[1] == 200),
              ('1 m, any height (front and left)', lambda d: where(d)[0] in ('front', 'left') and where(d)[1] == 100)]
    print('\nNumbers quoted in Sections 3.1 and 4: increase of H*(10) due to the floor (T/U - 1), range over the four '
          'radionuclides')
    for name, sel in groups:
        lo, hi = rng(sel)
        print(f'  {name:34s} {lo:5.1f} to {hi:5.1f} %')
    for n in C.NUCS:
        print(f'  {n}: 1 m in front {inc[n][LAB.index("front_g100_z110")]:.1f} %, 2 m in front '
              f'{inc[n][LAB.index("front_g200_z110")]:.1f} %, 1 m in front at 20 cm height '
              f'{inc[n][LAB.index("front_g100_z20")]:.1f} % and at 160 cm {inc[n][LAB.index("front_g100_z160")]:.1f} %')


if __name__ == '__main__':
    print('Supplementary runs R1 and R2 of article 1 - per GBq in the body, [177Lu]Lu-PSMA-617 distribution at 4 h')
    analyse_R1()
    analyse_R2()
    fn = os.path.join(_HERE, 'article1_supplementary_tests.csv')
    with open(fn, 'w', encoding='utf-8', newline='\n') as f:
        f.write('# R1: test,quantity,nuclide,detector,D,sD,G,sG,P,sP,D/G,s,G/P,s (per GBq in the body)\n')
        f.write('# R2: test,quantity,nuclide,detector,T,sT,F,sF,P,sP,T/U,s,T/P,s,U/P,s\n')
        f.write('\n'.join(CSV) + '\n')
    print(f'\nwritten: article1_supplementary_tests.csv ({len(CSV)} rows)')
