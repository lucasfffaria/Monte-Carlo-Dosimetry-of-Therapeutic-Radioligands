# -*- coding: utf-8 -*-
"""Article 1 - tables of the Supplementary Material (Supplementary_material_Article1_EANM.docx, Tables S1-S25).

Prints every table of the supplementary document in Markdown, followed by the numbers quoted in its text. All values
come from the files of this package:
  - model definitions: 3_model_scripts/data.py (emission data, energy grid, coefficients, regions), postproc_core.py
    (biokinetic model, region weights, folding), 1_MCNP_inputs/metadados_runs.json and massas_orgaos_fantomas.json;
  - results: article1_results.json and article1_detector_results.json (a1_03, a1_02), article1_cumulative_fraction.csv
    (a1_04), article1_supplementary_tests.csv (a1_05) and 2_MCNP_results/run_summary.csv;
  - quantities that no result file stores (ratios to ICRP 116; coefficients for other activity distributions; ratios of
    the simplified sources and contributions of the tissues to E; D/G by source region in run R1) are recomputed from
    the MCTAL files with the code of a1_01, a1_03 and a1_05, executed or imported here without printing and without
    writing any file (as a1_04 does with a1_03).
The script writes nothing; run it after a1_02-a1_05, whose output files it reads. The output is UTF-8.
"""
import csv
import io
import json
import os
import re
import runpy
import sys
from contextlib import redirect_stdout
from datetime import datetime

import numpy as np

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

# response matrices of a1_03 (same code, executed with printing suppressed, as in a1_04) and its time model
_SRC = open(os.path.join(_HERE, 'a1_03_article1_results.py'), encoding='utf-8').read()
_PRE, _POS = _SRC.split('# ------------------------------------------------------------------ time model')
exec(_PRE.replace("print(", "(lambda *a,**k: None)("))
exec('# ---- time model' + _POS.split("# H*(10) at 1 m and 2 m vs time")[0])
exec(re.search(r'^demir1=.*$', _POS, re.M).group(0))           # Demir et al. (2016), as printed by a1_03

import a1_05_supplementary_tests as T5                            # noqa: E402  (helpers of the supplementary runs)
from mctal import read_mctal                                      # noqa: E402

J = json.load(open(os.path.join(_HERE, 'article1_results.json'), encoding='utf-8'))
DET = json.load(open(os.path.join(_HERE, 'article1_detector_results.json'), encoding='utf-8'))
P0 = 'P0_paciente_em_pe'
DETS = {d[0]: d for d in META[P0]['detectors']}
NUC = {'Lu177': '¹⁷⁷Lu', 'Tb161': '¹⁶¹Tb', 'Pb212': '²¹²Pb', 'Ac225': '²²⁵Ac'}
EMIT = {'Lu-177': '¹⁷⁷Lu', 'Tb-161': '¹⁶¹Tb', 'Pb-212': '²¹²Pb', 'Bi-212': '²¹²Bi', 'Tl-208': '²⁰⁸Tl',
        'Ac-225': '²²⁵Ac', 'Fr-221': '²²¹Fr', 'Bi-213': '²¹³Bi', 'Tl-209': '²⁰⁹Tl'}
POS = (('FRONT', 'front', 'Face to face'), ('SIDE', 'left', 'Side by side'), ('BACK', 'back', 'Behind'))
GAPS = (30, 50, 100, 200)
TRUNK = np.array([0.0, 0.0, 105.0])                               # in-body point source (trunk centre)
QUOTED = []


def f3(x):
    """three significant figures; values from 1000 up as integers (as in the manuscript tables)"""
    if x == 0:
        return '0'
    a = abs(x)
    if a >= 1000:
        s = f'{x:.0f}'
    else:
        d = 2 - int(np.floor(np.log10(a)))
        s = f'{x:.{max(d, 0)}f}'
    return s.replace('-', '−')


def num(x, fmt):
    return format(x, fmt).replace('-', '−')


def pct(rel):
    return f'{100 * rel:.1f}'


def table(num_, title, header, rows):
    print(f'\n### Table {num_}. {title}\n')
    print('| ' + ' | '.join(header) + ' |')
    print('|' + '|'.join(['---'] * len(header)) + '|')
    for r in rows:
        assert len(r) == len(header), (num_, r)
        print('| ' + ' | '.join(str(x) for x in r) + ' |')


def ids_text(ids):
    """organ IDs as ranges: 1–8, 63, 65–71"""
    ids = sorted(ids)
    out, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(f'{ids[i]}' if i == j else f'{ids[i]}–{ids[j]}')
        i = j + 1
    return ', '.join(out)


def pretty(desc):
    return desc.replace(' & ', ' and ').replace('+', ' and ').replace('/', ' or ')


DIAG = {'45': 'left-back', '135': 'right-back', '225': 'right-front', '315': 'left-front'}


def label_text(lab):
    m = re.match(r'(front|back|left|right)_g(\d+)_z(\d+)', lab)
    if m:
        return f'{m.group(1)}, {m.group(2)} cm, z = {m.group(3)}'
    m = re.match(r'diag(\d+)_R(\d+)_z(\d+)', lab)
    return f'{DIAG[m.group(1)]} diagonal, R = {m.group(2)} cm, z = {m.group(3)}'


def det_order(lab):
    m = re.match(r'(front|back|left|right)_g(\d+)_z(\d+)', lab)
    if m:
        return (('front', 'back', 'left', 'right').index(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.match(r'diag(\d+)_R(\d+)_z(\d+)', lab)
    return (4 + (45, 135, 225, 315).index(int(m.group(1))), int(m.group(2)), 0)


def h_det(lab, n, w):
    return w @ RESP[('H', lab, n)][0] * GBQH


def e_pos(pos, gap, n, w):
    return w @ RESP[('Es', pos, gap, n)][0] * GBQH


print('Article 1 - Supplementary Material: Tables S1-S25 (Markdown)')
print('Rates per GBq in the body use the [177Lu]Lu-PSMA-617 distribution at 4 h, bladder just emptied, unless stated.')

# ================================================================================================ S1 phantoms
rows = []
for key, role in (('PAMstd', 'Patient (AM, standing)'), ('TAMstd', 'Exposed persons (AM)'),
                  ('TAFstd', 'Exposed persons (AF)')):
    m = MASS[key]
    vx = ' × '.join(f'{10 * v:.4g}' for v in m['vox_cm'])
    k, j, i = m['shape_kji']
    nel = sum(s['n_elements'] for s in m['slabs'])
    rows.append([key, role, vx, f'{i} × {j} × {k}', len(m['slabs']), f'{nel / 1e6:.2f}', f'{m["total_mass_g"] / 1e3:.1f}'])
table('S1', 'Voxel phantoms as implemented in MCNP6.3',
      ['Phantom', 'Role', 'Voxel size (mm)', 'Voxel array (columns × rows × slices)', 'Slabs',
       'Lattice elements (10⁶)', 'Mass (kg)'], rows)

# ================================================================================================ S2 regions
org = MASS['PAMstd']['organs']
rows = []
for r, (tag, desc, ids) in enumerate(D.REGIONS):
    if tag == 'SPO':
        ids_s, mass = 'spongiosa ' + ids_text(D.RBM_SITES) + ' (Table S3)', sum(org[str(s)]['mass_g'] for s in D.RBM_SITES)
    elif ids is None:
        ids_s, mass = '–', None
    else:
        ids_s, mass = ids_text(ids), sum(org[str(i)]['mass_g'] for i in ids if str(i) in org)
    rows.append([r + 1, tag, pretty(desc), ids_s, '–' if mass is None else f'{mass:.0f}',
                 f'{META[P0]["region_time_tag_shakes"][tag]:.0f}'])
table('S2', 'Source regions of the patient phantom',
      ['No.', 'Tag', 'Region', 'ICRP 110 organ IDs', 'Mass (g)', 'Birth-time tag (shakes)'], rows)

# ================================================================================================ S3 RBM sites
rows, tot = [], sum(D.SKEL_MASS[s]['AM_M'] for s in D.RBM_SITES)
for s in D.RBM_SITES:
    rows.append([org[str(s)]['name'], s, f'{org[str(s)]["mass_g"]:.0f}', f'{D.SKEL_MASS[s]["AM_M"]:.1f}',
                 f'{100 * D.SKEL_MASS[s]["AM_M"] / tot:.1f}'])
rows.append(['Total', '', f'{sum(org[str(s)]["mass_g"] for s in D.RBM_SITES):.0f}', f'{tot:.1f}', '100.0'])
table('S3', 'Spongiosa sites of the source region SPO and their weights (red-marrow mass of the adult male)',
      ['Bone site', 'ICRP 110 ID', 'Spongiosa mass in the patient (g)', 'Red-marrow mass, ICRP 116 (g)',
       'Share of the SPO activity (%)'], rows)

# ================================================================================================ S4 grid
emits = []
for nuc in NUCS:
    for name, mult, lines in D.NUCLIDE_DATA[nuc]['chain']:
        emits += [(EMIT[name], e) for e, y in lines if e >= 20.0]
rows = []
for k, e in enumerate(D.EGRID):
    co = sorted({f'{nm} {el:.1f}' for nm, el in emits if abs(el - 1000 * e) < 0.06})
    rows.append([k + 1, f'{1000 * e:.2f}', ', '.join(co) if co else '–', f'{100 * PE[k]:.2f}'] +
                [('–' if AN[n][k] == 0 else ('<0.0001' if AN[n][k] < 5e-5 else f'{AN[n][k]:.4f}')) for n in NUCS])
rows.append(['', 'Sum', '', '100.00'] + [f'{AN[n].sum():.4f}' for n in NUCS])
QUOTED.append(f'S4: {sum(1 for r in rows[:-1] if r[2] != "–")} of the {len(D.EGRID)} grid energies coincide with an emission line; '
              f'sampling probability {100 * PE.min():.2f}-{100 * PE.max():.2f} %; '
              f'{sum(len([1 for e, y in l if e >= 20.0]) for nuc in NUCS for _, _, l in D.NUCLIDE_DATA[nuc]["chain"])} emission lines in total')
table('S4', 'Source energy grid: sampling probability and photons per decay of the parent assigned to each energy',
      ['k', 'Eₖ (keV)', 'Coincident emission line (keV)', 'Sampling probability (%)'] + [NUC[n] for n in NUCS], rows)

# ================================================================================================ S5 emissions
rows = []
for nuc in NUCS:
    for name, mult, lines in D.NUCLIDE_DATA[nuc]['chain']:
        ll = sorted((e, y) for e, y in lines if e >= 20.0)
        txt = '; '.join(f'{e:.4f}'.rstrip('0').rstrip('.') + f' ({y:g})' for e, y in ll)
        rows.append([NUC[nuc], EMIT[name], f'{mult:g}', len(ll), txt])
table('S5', 'Photon emissions above 20 keV used for the four radionuclides',
      ['Radionuclide', 'Emitter', 'Decays per decay of the parent', 'Lines',
       'Energy in keV (photons per 100 decays of the emitter)'], rows)

# ================================================================================================ S6 activity
TIMES = (0.5, 4, 24, 72, 168)
W = {t: weights_from_model(lu_model(t, 0.0)) for t in TIMES}
rows = []
for r, (tag, desc, ids) in enumerate(D.REGIONS[:22]):
    rows.append([tag, pretty(desc.split(' (')[0])] + [f'{100 * W[t][r] / W[t].sum():.2f}' for t in TIMES])
rows.append(['Total', 'Activity in the body (% of the administered activity)'] + [f'{lu_model(t)["WB"]:.1f}' for t in TIMES])
table('S6', 'Distribution of the activity in the body over the source regions (% of the activity in the body), '
      'bladder just emptied', ['Tag', 'Region'] + [f'{t:g} h' for t in TIMES], rows)

# ================================================================================================ S7 time dependence
W05, W72 = (W[t] / W[t].sum() for t in (0.5, 72))
rows = []
for n in NUCS:
    def worst(items):
        lab, a, b = max(items, key=lambda x: abs(x[2] / x[1] - 1))
        return f'{num(100 * (b / a - 1), "+.1f")} ({lab})'
    fr = [(f'{g} cm', h_det(f'front_g{g}_z110', n, W05), h_det(f'front_g{g}_z110', n, W72)) for g in (10, 30, 50, 100, 200, 300)]
    al = [(f'{f} {g} cm', h_det(f'{f}_g{g}_z110', n, W05), h_det(f'{f}_g{g}_z110', n, W72))
          for f in ('front', 'back', 'left', 'right') for g in (10, 30, 50, 100, 200, 300)]
    ee = [(f'{nm.lower()} {g} cm', e_pos(p, g, n, W05), e_pos(p, g, n, W72)) for p, _, nm in POS for g in GAPS]
    fb = [x for x in ee if not x[0].startswith('side')]
    rows.append([NUC[n], worst(fr), worst(al), worst(fb), worst(ee)])
    QUOTED.append(f'S7: {n}, change from the 0.5 h to the 72 h distribution, largest |change|: H*(10) front '
                  f'{max(abs(b / a - 1) for _, a, b in fr) * 100:.1f} %, H*(10) all 24 detectors at 110 cm '
                  f'{max(abs(b / a - 1) for _, a, b in al) * 100:.1f} %, E face to face and behind '
                  f'{max(abs(b / a - 1) for _, a, b in fb) * 100:.1f} %, E all positions {max(abs(b / a - 1) for _, a, b in ee) * 100:.1f} %')
table('S7', 'Largest change of the coefficients per GBq in the body when the 72 h distribution replaces the 0.5 h '
      'distribution (%, where it occurs)',
      ['Radionuclide', 'H*(10) in front, 10 cm to 3 m', 'H*(10) at the 24 detectors at 110 cm height',
       'E, face to face and behind', 'E, all positions'], rows)

# ================================================================================================ S8 coefficients
hs = dict(zip(D.HSTAR_E, D.HSTAR_F))
ka = dict(zip(D.KA_E, D.KA_F))
rows = [[f'{e:g}', f'{hs[e]:g}' if e in hs else '–', f'{ka[e]:.5g}' if e in ka else '–'] for e in sorted(set(hs) | set(ka))]
table('S8', 'Fluence-to-dose coefficients entered on the DE/DF cards (log-log interpolation)',
      ['Photon energy (MeV)', 'H*(10)/Φ, ICRP 74 (pSv cm²)', 'Ka/Φ, derived from ICRP 116 (pGy cm²)'], rows)

# ================================================================================================ S9 tissues
names = {t[0]: (t[1], t[2], t[3]) for t in D.TISSUES}
mam, maf = MASS['TAMstd']['tissues'], MASS['TAFstd']['tissues']
order = ['GON', 'RBM', 'COL', 'LUN', 'STW', 'BRE', 'UBW', 'OES', 'LIV', 'THY', 'BSE', 'BRA', 'SAL', 'SKI'] + REM + ['LEN', 'REC']
rows = []
for t in order:
    if t in ('RBM', 'BSE'):                      # spongiosa approximation: mass-weighted dose to the spongiosa sites
        kk = ('AM_M', 'AM_F') if t == 'RBM' else ('TM50_M', 'TM50_F')
        sites = [s for s in D.SPONGIOSA if s in D.SKEL_MASS]
        desc = 'Red bone marrow' if t == 'RBM' else 'Endosteum (bone surface)'
        rows.append([desc, f'{WT[t]:g}', 'spongiosa ' + ids_text(sites), 'same',
                     f'{sum(D.SKEL_MASS[s][kk[0]] for s in sites):.0f}', f'{sum(D.SKEL_MASS[s][kk[1]] for s in sites):.0f}'])
        continue
    desc, am_ids, af_ids = names[t]
    w = f'{WT[t]:g}' if t in WT else (f'{WT["REM"]:g}/13' if t in REM else '–')
    desc += '' if t in WT else (' (remainder)' if t in REM else ' (not in E)')
    rows.append([desc, w, ids_text(am_ids), 'same' if af_ids == am_ids else ids_text(af_ids),
                 f'{mam[t]["mass_g"]:.1f}', f'{maf[t]["mass_g"]:.1f}'])
table('S9', 'Target tissues of the exposed phantoms: ICRP 103 tissue weighting factors, ICRP 110 organ IDs and masses',
      ['Tissue', 'w_T', 'Organ IDs, AM', 'Organ IDs, AF', 'Mass, AM (g)', 'Mass, AF (g)'], rows)

# ================================================================================================ S10 detectors
groups = {}
for lab, xyz, rad in META[P0]['detectors']:
    m = re.match(r'(front|back|left|right)_g(\d+)_z(\d+)', lab)
    key = (m.group(1), m.group(2)) if m else re.match(r'(diag\d+)_R(\d+)', lab).groups()
    groups.setdefault(key, []).append((xyz, rad))
rows = []
for key in sorted(groups, key=lambda k: det_order(f'{k[0]}_g{k[1]}_z0' if not k[0].startswith('diag')
                                                  else f'{k[0]}_R{k[1]}_z0')):
    pts = groups[key]
    x, y = pts[0][0][0], pts[0][0][1]
    d110 = np.linalg.norm(np.array([x, y, 110.0]) - TRUNK)
    if key[0].startswith('diag'):
        lab, gap = f'{DIAG[key[0][4:]].capitalize()} diagonal ({key[0][4:]}°)', f'R = {key[1]}'
    else:
        lab, gap = key[0].capitalize(), key[1]
    rows.append([lab, gap, num(x, '.1f'), num(y, '.1f'), ', '.join(f'{z:g}' for z in sorted(p[0][2] for p in pts)),
                 f'{pts[0][1]:g}', f'{d110:.1f}'])
env = [round(v, 1) + 0.0 for v in META[P0]['placements'][0]['envelope'][0]]
QUOTED.append(f'S10: bounding box of the patient: x {env[0]:.1f} to {env[1]:.1f} cm, y {env[2]:.1f} to {env[3]:.1f} cm, '
              f'z {env[4]:.1f} to {env[5]:.1f} cm; trunk centre {TRUNK.tolist()}')
table('S10', 'Air-sphere detectors of the patient-alone run (76 spheres)',
      ['Direction', 'Gap from the bounding box (cm)', 'x (cm)', 'y (cm)', 'Heights z (cm)', 'Sphere radius (cm)',
       'Distance from the trunk centre at z = 110 cm (cm)'], rows)

# ================================================================================================ S11 ICRP 116
_buf = io.StringIO()
with redirect_stdout(_buf):
    V1 = runpy.run_path(os.path.join(_HERE, 'a1_01_verification_icrp116.py'))
rows = []
for i, en in enumerate(V1['En']):
    rows.append([f'{1000 * en:g}'] + [f'{V1["R"][g][1][i]:.3f} ({V1["R"][g][0][i]:.3f})' for g in ('AP', 'PA', 'LLAT', 'RLAT')])
table('S11', 'Effective dose per fluence relative to ICRP 116: spongiosa approximation (dose-response functions)',
      ['Energy (keV)', 'AP', 'PA', 'LLAT', 'RLAT'], rows)
sel = np.array(V1['En']) >= 0.03
alls = np.array([V1['R'][g][1] for g in ('AP', 'PA', 'LLAT', 'RLAT')])
alld = np.array([V1['R'][g][0] for g in ('AP', 'PA', 'LLAT', 'RLAT')])
QUOTED.append(f'S11: spongiosa method, 30 keV-3 MeV: mean {alls[:, sel].mean():.3f}, range {alls[:, sel].min():.3f}-'
              f'{alls[:, sel].max():.3f}; 20 keV: {alls[:, ~sel].min():.3f}-{alls[:, ~sel].max():.3f}; DRF method below '
              f'100 keV: {alld[:, np.array(V1["En"]) < 0.1].min():.3f}-{alld[:, np.array(V1["En"]) < 0.1].max():.3f}')

# ================================================================================================ S12 time curve
rows = []
for tc in J['timecurve']:
    t = tc['t']
    d1, d2 = (f'{demir1[0]} (start)', f'{demir2[0]} (start)') if t < 1 else (demir1.get(int(t), '–'), demir2.get(int(t), '–'))
    rows.append([f'{t:g}', f'{lu_model(t)["WB"]:.1f}', f3(tc['H05_empty']), f'{f3(tc["H1_empty"])}–{f3(tc["H1_full"])}',
                 f'{f3(tc["H2_empty"])}–{f3(tc["H2_full"])}', f'{d1}', f'{d2}'])
table('S12', 'Modelled H*(10) rate (µSv h⁻¹) in front of a patient given 7.4 GBq [¹⁷⁷Lu]Lu-PSMA-617, 110 cm above the '
      'floor, and the means measured by Demir et al. (2016)',
      ['Time after administration (h)', 'Activity in the body (% IA)', '0.5 m, bladder empty', '1 m, bladder empty–full',
       '2 m, bladder empty–full', 'Demir et al., 1 m', 'Demir et al., 2 m'], rows)

# ================================================================================================ S13-S15 R1 and R2
S = {}
with open(os.path.join(_HERE, 'article1_supplementary_tests.csv'), encoding='utf-8') as f:
    for line in f:
        if not line.startswith('#'):
            p = line.strip().split(',')
            S[(p[0], p[1], p[2], p[3])] = [float(x) for x in p[4:]]
rows = []
for lab in T5.KEY:
    h, k = S[('R1', 'H*(10)', 'Tb161', lab)], S[('R1', 'Ka', 'Tb161', lab)]
    rows.append([label_text(lab), f'{h[6]:.4f} ± {h[7]:.4f}', f'{k[6]:.4f} ± {k[7]:.4f}', f'{h[8]:.3f} ± {h[9]:.3f}'])
table('S13', 'Run R1, ¹⁶¹Tb: ratio of the result with the real emission lines (D) to that with the lines folded on '
      'the energy grid (G), and ratio of G to the production run P0 (P)',
      ['Detector (direction, gap, height in cm)', 'D/G, H*(10)', 'D/G, Ka', 'G/P, H*(10)'], rows)
# D/G by source region for all photons leaving the patient (tally F344), as in a1_05
T1 = read_mctal(T5.C.paths.mctal(T5.R1))['tallies']
m1 = T5.META[T5.R1]
ent = m1['energy_entries']
pp = np.array([e['p_deck'] for e in ent])
pp = pp / pp.sum()
fac = np.array([e['photons_per_decay'] for e in ent]) / pp
isl = np.array([e['kind'] == 'line' for e in ent])
v344, e344 = T5.cube(T1[344])
nE, nR = len(ent), len(m1['regions'])


def dg(it):
    v = v344[0, 0, :nE, 0, 0, 0, 0, it] * fac
    s = v * e344[0, 0, :nE, 0, 0, 0, 0, it]
    return T5.ratio(v[isl].sum(), np.sqrt((s[isl] ** 2).sum()), v[~isl].sum(), np.sqrt((s[~isl] ** 2).sum()))


cells = []
for k, reg in enumerate(m1['regions']):
    r, sr = dg(k)
    if np.isfinite(r) and r > 0:
        cells.append((reg[0], f'{r:.4f} ± {sr:.4f}'))
rt, srt = dg(nR)
cells.append(('All regions', f'{rt:.5f} ± {srt:.5f}'))
while len(cells) % 3:
    cells.append(('', ''))
rows = [[cells[i][0], cells[i][1], cells[i + 1][0], cells[i + 1][1], cells[i + 2][0], cells[i + 2][1]]
        for i in range(0, len(cells), 3)]
table('S14', 'Run R1: D/G for all photons leaving the patient (tally F344), by source region',
      ['Region', 'D/G', 'Region', 'D/G', 'Region', 'D/G'], rows)
QUOTED.append(f'S14: all photons leaving the patient, D/G = {rt:.5f} +- {srt:.5f}')
R2KEY = T5.KEY + [f'front_g100_z{z}' for z in T5.HEIGHTS if z != 110] + [f'left_g100_z{z}' for z in T5.HEIGHTS if z != 110]
rows = []
for lab in R2KEY:
    rows.append([label_text(lab)] + [f'{100 * (S[("R2", "H*(10)", n, lab)][6] - 1):.1f} ± {100 * S[("R2", "H*(10)", n, lab)][7]:.1f}'
                                     for n in NUCS])
table('S15', 'Run R2: increase of H*(10) due to a 30 cm concrete floor, T/U − 1 (%), with its statistical uncertainty',
      ['Detector (direction, gap, height in cm)'] + [NUC[n] for n in NUCS], rows)
for z in (90, 110, 135):
    a, b = S[('R2', 'H*(10)', 'Lu177', f'front_g200_z{z}')], S[('R2', 'H*(10)', 'Lu177', f'front_g100_z{z}')]
    QUOTED.append(f'S15: 177Lu, z = {z} cm, H*(10)(2 m)/H*(10)(1 m) in front: without floor (P0) {a[4] / b[4]:.3f}, '
                  f'with floor (T) {a[0] / b[0]:.3f}, without the floor-scattered part (U) {(a[0] - a[2]) / (b[0] - b[2]):.3f}')
up = {n: [S[('R2', 'H*(10)', n, lab)][10] for lab in T5.LAB] for n in NUCS}
QUOTED.append('S15: U/P over the 76 detectors (mean): ' + '; '.join(f'{n} {min(up[n]):.3f}-{max(up[n]):.3f} ({np.mean(up[n]):.3f})' for n in NUCS))

# ================================================================================================ S16 H*(10)
labs = sorted(DET, key=det_order)
rows = [[label_text(lab)] + [f'{f3(DET[lab][n]["H"])} ({pct(DET[lab][n]["rel"])})' for n in NUCS] for lab in labs]
table('S16', 'H*(10) rate per GBq in the body (µSv h⁻¹ GBq⁻¹) at the 76 detectors of the patient-alone run '
      '(relative statistical uncertainty, %)', ['Detector (direction, gap, height in cm)'] + [NUC[n] for n in NUCS], rows)
QUOTED.append('S16: relative uncertainty of H*(10) over the 76 detectors: ' +
              '; '.join(f'{n} {100 * min(DET[l][n]["rel"] for l in labs):.1f}-{100 * max(DET[l][n]["rel"] for l in labs):.1f} %'
                        for n in NUCS))

# ================================================================================================ S17 air kerma
rows = []
for n in NUCS:
    rat = np.array([DET[lab][n]['H'] / DET[lab][n]['K'] for lab in labs])
    r1 = DET['front_g100_z110'][n]['H'] / DET['front_g100_z110'][n]['K']
    rows.append([NUC[n]] + [f3(DET[f'{f}_g100_z110'][n]['K']) for f in ('front', 'back', 'left', 'right')] +
                [f'{r1:.3f}', f'{rat.min():.3f}–{rat.max():.3f}'])
    QUOTED.append(f'S17: {n} H*(10)/Ka at the 76 detectors relative to the value 1 m in front: '
                  f'{100 * (rat.min() / r1 - 1):+.1f} % to {100 * (rat.max() / r1 - 1):+.1f} %; detectors beyond 3 %: '
                  + (', '.join(l for l in labs if abs(DET[l][n]['H'] / DET[l][n]['K'] / r1 - 1) > 0.03) or 'none'))
table('S17', 'Air-kerma rate per GBq in the body 1 m from the patient at 110 cm height (µGy h⁻¹ GBq⁻¹) and ratio '
      'H*(10)/Ka (Sv Gy⁻¹)',
      ['Radionuclide', 'Ka, front', 'Ka, back', 'Ka, left', 'Ka, right', 'H*(10)/Ka, 1 m in front',
       'H*(10)/Ka, range over the 76 detectors'], rows)

# ================================================================================================ S18 simplified sources (H*)
SEL = ['front_g10_z110', 'front_g30_z110', 'front_g50_z110', 'front_g100_z110', 'front_g200_z110', 'front_g300_z110',
       'back_g10_z110', 'back_g100_z110', 'back_g300_z110', 'left_g10_z110', 'left_g100_z110', 'left_g300_z110',
       'front_g100_z45', 'front_g100_z160']
rows = []
for lab in SEL:
    r = np.linalg.norm(np.array(DETS[lab][1]) - TRUNK)
    cells = []
    for n in NUCS:
        d = DET[lab][n]
        bare = J['Gamma'][n]['Hstar'] * (100.0 / r) ** 2
        cells.append(f'{d["Hp"] / d["H"]:.2f} / {d["Hl"] / d["H"]:.2f} / {bare / d["H"]:.2f}')
    rows.append([label_text(lab), f'{r:.1f}'] + cells)
table('S18', 'H*(10) of the simplified sources relative to the voxel source: in-body point / in-body line / bare point source',
      ['Detector (direction, gap, height in cm)', 'Distance from the trunk centre (cm)'] + [NUC[n] for n in NUCS], rows)

# ================================================================================================ S19 E
rows = []
for pos, face, name in POS:
    for g in GAPS:
        rows.append([name, f'{g}'] + [f'{f3(J["E"][f"{pos}_{g}_{n}"]["spong"])} ({pct(J["E"][f"{pos}_{g}_{n}"]["relerr"])})'
                                      for n in NUCS])
table('S19', 'Effective dose rate per GBq in the body (µSv h⁻¹ GBq⁻¹) of the exposed persons, sex-averaged, spongiosa '
      'approximation (conservative relative statistical uncertainty, %)', ['Position', 'Gap (cm)'] + [NUC[n] for n in NUCS], rows)
for n in NUCS:
    dr = [J['E'][f'{p}_{g}_{n}']['drf'] / J['E'][f'{p}_{g}_{n}']['spong'] for p, _, _ in POS for g in GAPS]
    QUOTED.append(f'S19: {n} E with dose-response functions / spongiosa: {min(dr):.3f}-{max(dr):.3f}')
dev = max(abs(e_pos(p, g, n, w4n) / J['E'][f'{p}_{g}_{n}']['spong'] - 1) for p, _, _ in POS for g in GAPS for n in NUCS)
QUOTED.append(f'check: E recomputed here / article1_results.json - 1: max {dev:.1e}')

# ================================================================================================ S20 E simplified
wp, wl = w_point(), w_line()
rows = []
for pos, face, name in POS:
    for g in GAPS:
        rows.append([name, f'{g}'] + [f'{(wp @ RESP[("Es", pos, g, n)][0]) / (w4n @ RESP[("Es", pos, g, n)][0]):.2f} / '
                                      f'{(wl @ RESP[("Es", pos, g, n)][0]) / (w4n @ RESP[("Es", pos, g, n)][0]):.2f}' for n in NUCS])
table('S20', 'Effective dose from the in-body point and line sources relative to the voxel source (point / line)',
      ['Position', 'Gap (cm)'] + [NUC[n] for n in NUCS], rows)

# ================================================================================================ S21 tissues
TNAME = {'GON': 'Gonads', 'RBM': 'Red bone marrow', 'COL': 'Colon', 'LUN': 'Lungs', 'STW': 'Stomach', 'BRE': 'Breast',
         'REM': 'Remainder (13 tissues)', 'UBW': 'Urinary bladder', 'OES': 'Oesophagus', 'LIV': 'Liver', 'THY': 'Thyroid',
         'BSE': 'Endosteum', 'BRA': 'Brain', 'SAL': 'Salivary glands', 'SKI': 'Skin'}
am100, af100 = target_organs_s('T03_AM_g100', 'TAMstd'), target_organs_s('T07_AF_g100', 'TAFstd')
contrib = {}
for pos in ('FRONT', 'SIDE', 'BACK'):
    dm, df = dict(am100[pos + 'AM']), dict(af100[pos + 'AF'])
    for d in (dm, df):
        d['RBM'], d['BSE'] = d['RBMs'], d['BSEs']
    for n in NUCS:
        a = AN[n]
        ht = {t: w4n @ ((0.5 * (dm[t][0] + df[t][0])) @ a) * GBQH for t in WT if t != 'REM'}
        ht['REM'] = np.mean([w4n @ ((0.5 * (dm[x][0] + df[x][0])) @ a) * GBQH for x in REM])
        e_tot = sum(WT[t] * ht[t] for t in WT)
        assert abs(e_tot / J['E'][f'{pos}_100_{n}']['spong'] - 1) < 1e-9
        contrib[(pos, n)] = {t: 100 * WT[t] * ht[t] / e_tot for t in WT}
rows = [[TNAME[t], f'{WT[t]:g}'] + [f'{contrib[(p, n)][t]:.1f}' for p in ('FRONT', 'SIDE') for n in NUCS]
        for t in sorted(WT, key=lambda t: -contrib[('FRONT', 'Lu177')][t])]
table('S21', 'Contribution of the tissues to the effective dose of the exposed person at a 1 m gap (% of E)',
      ['Tissue', 'w_T'] + [f'Face to face, {NUC[n]}' for n in NUCS] + [f'Side by side, {NUC[n]}' for n in NUCS], rows)
for p in ('FRONT', 'SIDE', 'BACK'):
    QUOTED.append(f'S21: {p} 1 m, breast + skin + thyroid + salivary glands: ' +
                  ', '.join(f'{n} {sum(contrib[(p, n)][t] for t in ("SKI", "BRE", "THY", "SAL")):.1f} %' for n in NUCS) +
                  '; gonads + stomach + urinary bladder: ' +
                  ', '.join(f'{n} {sum(contrib[(p, n)][t] for t in ("GON", "STW", "UBW")):.1f} %' for n in NUCS))

# ================================================================================================ S22 H*/E
rows = []
for pos, face, name in POS:
    for g in GAPS:
        lab = f'{face}_g{g}_z110'
        r = np.linalg.norm(np.array(DETS[lab][1]) - TRUNK)
        rows.append([name, f'{g}'] + [f'{h_det(lab, n, w4n) / e_pos(pos, g, n, w4n):.2f}' for n in NUCS] +
                    [f'{J["Gamma"][n]["Hstar"] * (100.0 / r) ** 2 / e_pos(pos, g, n, w4n):.1f}' for n in NUCS])
table('S22', 'Ratio of H*(10) to the effective dose of the exposed person: H*(10) at the survey-meter point of the '
      'patient-alone run, and H*(10) of the bare point source at the same point',
      ['Position', 'Gap (cm)'] + [f'Survey meter, {NUC[n]}' for n in NUCS] + [f'Bare point, {NUC[n]}' for n in NUCS], rows)

# ================================================================================================ S23 cumulative
cum = {n: [] for n in NUCS}
tt = []
with open(os.path.join(_HERE, 'article1_cumulative_fraction.csv'), encoding='utf-8') as f:
    for row in csv.DictReader(f):
        tt.append(float(row['t_after_release_h']))
        for n in NUCS:
            cum[n].append(float(row[n]))
TS = (2, 6, 12, 24, 48, 72, 168, 336)
rows = [[NUC[n]] + [f'{100 * np.interp(t, tt, cum[n]):.1f}' for t in TS] for n in NUCS]
table('S23', 'Cumulative share of the 30-day effective dose 1 m in front of the patient (%), release 4 h after '
      'administration, voiding every 4 h', ['Radionuclide'] + [f'{t:g} h' if t < 24 else f'{t // 24:g} d' for t in TS], rows)

# ================================================================================================ S24 release
rows = []
for n in NUCS:
    a0 = J['integrals'][f'{n}_4']['A0']
    h2 = [h_det(f'{f}_g200_z110', n, w4n) for f in ('front', 'back', 'left', 'right')]
    amax = J['release2m'][n]['Amax']
    h2rel = rate(('H', 'front_g200_z110', n), n, 4.0, 4.0) * a0          # at release, typical activity
    einf = J['integrals'][f'{n}_4']['E1m_front']
    rows.append([NUC[n]] + [f3(x) for x in h2] + [f'{a0:g}', f3(amax), f'{amax / a0:.3g}', f3(h2rel), f3(einf),
                                                 f'{einf / h2rel:.0f}'])
table('S24', 'Release criterion of 30 µSv h⁻¹ at 2 m: H*(10) at 2 m per GBq in the body, activity in the body that '
      'meets the criterion, and time-integrated E per unit H*(10) rate at release (4 h)',
      ['Radionuclide', 'H*(10) 2 m front (µSv h⁻¹ GBq⁻¹)', '2 m back', '2 m left', '2 m right', 'Typical activity (GBq)',
       'Activity at 30 µSv h⁻¹ (GBq)', 'Ratio to typical activity', 'H*(10) 2 m front at release (µSv h⁻¹)',
       'E, 1 m front, release to 30 d (µSv)', 'E per unit H*(10) rate (µSv per µSv h⁻¹)'], rows)

# ================================================================================================ S25 runs
DESC = {'V': 'Verification, broad parallel beam', 'P': 'Patient alone', 'T': 'Patient and three exposed persons',
        'R': 'Supplementary run'}
with open(os.path.join(T5.C.paths.RESULTS, 'run_summary.csv'), encoding='utf-8') as f:
    runs = list(csv.DictReader(f))
rows = []
for r in sorted(runs, key=lambda r: ('VPTR'.index(r['run'][0]), r['run'])):
    t0, t1 = (datetime.strptime(r[k], '%m/%d/%y %H:%M:%S') for k in ('start', 'end'))
    rows.append([r['run'], DESC[r['run'][0]], f'{int(r["nps"]):,}'.replace(',', '\u00a0'), f'{float(r["computer_time_min"]):.0f}',
                 f'{(t1 - t0).total_seconds() / 60:.0f}', r['warnings'], r['lost_particles'],
                 f'{r["tallies_passing_10_checks"]} of {r["tallies_total"]}'])
table('S25', 'Production and supplementary runs: histories, computer time, warnings, lost particles and statistical checks',
      ['Run', 'Configuration', 'Histories', 'Computer time (min)', 'Wall-clock time (min)', 'Warnings', 'Lost particles',
       'Tallies passing all ten checks'], rows)

print('\nNumbers quoted in the text of the Supplementary Material')
for q in QUOTED:
    print('  ' + q)
