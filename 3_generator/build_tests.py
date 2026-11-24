# -*- coding: utf-8 -*-
"""Decks of the two supplementary runs of article 1 (Sections 2.7, 3.1 and 4; Supplementary Sections S3.3 and S3.4), built with
the same modules as the production decks.

R1_P0_Tb161_direto.i
    Geometry and detectors of P0 (patient alone in air). Source regions sampled with the reference
    [177Lu]Lu-PSMA-617 activity distribution at 4 h, just after voiding (the distribution of the per-GBq results), and
    source energies in two groups sampled in the same run: (a) the 21 photon lines of 161Tb (half of the histories);
    (b) the same lines folded on the 27-energy grid with the linear-in-log-E weights of the post-processing (the other
    half). Tallies are binned by energy entry (FT SCX 6) and by region (T0), and an extra tally (F344) scores all
    photons that leave the patient. The ratio (a)/(b) measures the folding error.
R2_P0_piso_concreto.i
    P0 with a 30 cm floor of NIST Portland concrete (2.30 g/cm3) just below the soles. Source and tallies identical to
    those of P0 (any radionuclide and activity distribution by the same folding) plus CF cards: the second (flagged)
    print of each tally is the contribution of the photons that have left the floor.

Usage: python build_tests.py [output folder]     (default: work/decks, as for build_all.py)
Writes the two decks and metadados_testes.json (source energies and regions of R1, floor of R2), first checking that
P0 is rebuilt identical byte for byte to the official deck, and then runs the deck checks of validate.py.
On 4 November 2026 the output was identical byte for byte to the files in 1_MCNP_inputs. The decks read the lattice
file lat_PAMstd.txt, so they are run from 1_MCNP_inputs (07_TESTES_R1_R2.bat).
"""
import json
import os
import sys
import tempfile

import numpy as np

import paths
import data as D
import scenes
from mcnpgen import (Deck, wrap, fsig, AIR_MAT, AIR_RHO, EDIST, ERG0,
                     patient_source, detector_tallies, run_cards)
sys.path.insert(0, str(paths.ANALYSIS))
import postproc_core as C         # noqa: E402  (biokinetic model, region weights, line folding)

NPS = 2_000_000_000
CTME = 1200          # minutes of computer time summed over the threads: about 1 h of wall-clock time on 20 threads
NUC = 'Tb161'
F_LINES = 0.5        # fraction of the histories that sample the real lines (the others sample the grid)
FLOOR_CELL, FLOOR_MAT, FLOOR_RHO = 7900, 300, 2.30
FLOOR_BOX = (-420.0, 420.0, -420.0, 420.0, -30.01, -0.01)     # cm; top 0.1 mm below the soles (z = 0)
CONCRETE = [(1000, 0.010000), (6000, 0.001000), (8000, 0.529107), (11000, 0.016000), (12000, 0.002000),
            (13000, 0.033872), (14000, 0.337021), (19000, 0.013000), (20000, 0.044000), (26000, 0.014000)]
SPECTRUM_DETECTORS = ('front_g100_z110', 'left_g100_z110', 'back_g100_z110')
DECKS = {'R1': 'R1_P0_Tb161_direto', 'R2': 'R2_P0_piso_concreto'}
META_FILE = 'metadados_testes.json'


class FloorDeck(Deck):
    """P0 deck with the floor material appended to the material list"""
    def materials(self):
        out = super().materials()
        toks = []
        for z, f in CONCRETE:
            toks += [str(z), f'-{f:.6f}']
        out.append('c ---- floor: NIST Portland concrete, 2.30 g/cm3')
        return out + wrap(toks, first=f'm{FLOOR_MAT}')


def p0_scene(deck):
    """patient and detectors exactly as in scenes.deck_P0"""
    P = scenes.patient_standing(deck)
    cells, labels = [], []
    for n, (lbl, c, r) in enumerate(scenes.p0_detectors()):
        s = deck.sph(c, r, lbl)
        cn = 7101 + n
        deck.det_vol[cn] = 4 / 3 * np.pi * r ** 3
        deck.cell(cn, AIR_MAT, AIR_RHO, [f'-{s}'], comment=lbl, outside_tokens=[str(s)])
        cells.append(cn)
        labels.append(lbl)
    return P, cells, labels


def region_weights_4h():
    """fraction of the activity in the body per region at 4 h, just after voiding (normalised to 1)"""
    w = C.weights_from_model(C.lu_model(4.0, 0.0))
    return w / w.sum()


def energy_entries():
    """energy entries of R1: the 21 lines of 161Tb and the same lines folded on the grid"""
    lines = D.nuclide_lines(NUC)
    Y = sum(y for _, y in lines)
    a = C.line_matrix(NUC)
    assert abs(a.sum() - Y) < 1e-12
    ent = [dict(kind='line', E_MeV=e, photons_per_decay=y, p=F_LINES * y / Y) for e, y in lines]
    ent += [dict(kind='grid', E_MeV=float(D.EGRID[k]), grid_index=k, photons_per_decay=float(a[k]),
                 p=(1 - F_LINES) * float(a[k]) / Y) for k in range(len(a)) if a[k] > 0]
    for e in ent:                     # probability as written in the deck (6 significant digits)
        e['p_deck'] = float(fsig(e['p'], 6))
    return ent, Y


def replace_energy_cards(deck, ent):
    i0 = next(i for i, l in enumerate(deck.data) if l.startswith('c ---- source energy grid'))
    i1 = next(i for i, l in enumerate(deck.data) if l.startswith('c ---- source region tags'))
    nl = sum(e['kind'] == 'line' for e in ent)
    new = [f'c ---- source energies (MeV): entries 1-{nl} = the {nl} 161Tb photon lines',
           f'c      (half of the histories); entries {nl + 1}-{len(ent)} = the same lines folded on',
           'c      the energy grid with linear-in-log-E weights (other half). Tallies are',
           f'c      binned by entry with FT SCX {EDIST}.']
    new += wrap([fsig(e['E_MeV'], 6) for e in ent], first=f'SI{EDIST} L')
    new += wrap([fsig(e['p_deck'], 6) for e in ent], first=f'SP{EDIST}')
    deck.data[i0:i1] = new


def world_tally():
    L = ['c ---- H*(10)-weighted track length in the world air cell: all photons that',
         'c      leave the patient (angle-integrated check of the energy folding)',
         'F344:P 9998', 'SD344 1', 'FC344 H*(10)-weighted track length in the world air cell (rel. units)']
    L += wrap([fsig(e, 4) for e in D.HSTAR_E], first='DE344 LOG')
    L += wrap([fsig(v, 4) for v in D.HSTAR_F], first='DF344 LOG')
    L.append(f'FT344 SCX {EDIST}')
    return L


def source_meta(d):
    ent = d.src_entries
    tot = sum(e['p'] for e in ent)
    return dict(regions=[(r[0], r[1]) for r in d.regions],
                region_prob={r[0]: sum(e['p'] for e in ent if e['ed'] == ERG0 + i + 1) / tot
                             for i, r in enumerate(d.regions)},
                region_time_tag_shakes={r[0]: 1000.0 * (i + 1) for i, r in enumerate(d.regions)},
                n_source_entries=len(ent), source_stats=d.src_stats)


def deck_R1():
    d = Deck('R1: P0 geometry, 161Tb photon lines versus their grid folding (test)', world_radius=600.0)
    P, cells, labels = p0_scene(d)
    d.finish_world()
    w = region_weights_4h()
    regs = [r for r, x in zip(D.REGIONS, w) if x > 0]
    patient_source(d, P, region_weights=[float(x) for x in w if x > 0], line_z=(2.0, 174.0), regions=regs)
    ent, Y = energy_entries()
    replace_energy_cards(d, ent)
    d.notes = ['Review test R1: patient alone in air (geometry and detectors of P0).',
               'Source regions sampled with the reference [177Lu]Lu-PSMA-617 activity distribution at 4 h '
               '(just after voiding; no bladder contents, no point or line source).',
               'Energies: half of the histories sample the 21 photon lines of 161Tb, the other half the same '
               'lines folded on the 27-energy grid of the production runs. Tallies binned by source energy entry '
               '(FT SCX 6) and by region (T0); their ratio measures the folding error.']
    spec = [cells[labels.index(x)] for x in SPECTRUM_DETECTORS]
    detector_tallies(d, cells, labels, spectrum_cells=spec)
    d.data += world_tally()
    run_cards(d, NPS, CTME)
    m = dict(kind='R1', nuclide=NUC, f_lines=F_LINES, photons_per_decay_total=Y, energy_entries=ent,
             region_weights={r[0]: float(x) for r, x in zip(D.REGIONS, w) if x > 0},
             tally_index={'F304': labels, 'F314': labels, 'F324': list(SPECTRUM_DETECTORS), 'F344': ['world_air']},
             **source_meta(d))
    return d, m


def deck_R2():
    d = FloorDeck('R2: P0 with a 30 cm concrete floor under the patient (test)', world_radius=600.0)
    P, cells, labels = p0_scene(d)
    s = d.rpp(FLOOR_BOX, 'concrete floor')
    d.cell(FLOOR_CELL, FLOOR_MAT, FLOOR_RHO, [f'-{s}'], comment='concrete floor', outside_tokens=[str(s)])
    d.finish_world()
    patient_source(d, P, line_z=(2.0, 174.0))
    d.notes = ['Review test R2: P0 (patient alone, same source and tallies) plus a 30 cm floor of NIST '
               'Portland concrete (2.30 g/cm3) from z = -30.01 to -0.01 cm under the soles (z = 0), '
               f'{FLOOR_BOX[1] - FLOOR_BOX[0]:.0f} x {FLOOR_BOX[3] - FLOOR_BOX[2]:.0f} cm.',
               'CF cards: the second (flagged) print of each tally is the contribution of photons that have '
               'left the floor cell, i.e. the floor-scattered radiation.'] + d.notes
    spec = [cells[labels.index(x)] for x in SPECTRUM_DETECTORS]
    detector_tallies(d, cells, labels, spectrum_cells=spec)
    d.data += ['c ---- cell flagging: contributions of photons that have left the floor',
               f'CF304 {FLOOR_CELL}', f'CF314 {FLOOR_CELL}', f'CF324 {FLOOR_CELL}']
    run_cards(d, NPS, CTME)
    m = dict(kind='R2', floor=dict(cell=FLOOR_CELL, material=FLOOR_MAT, density_g_cm3=FLOOR_RHO, box_cm=FLOOR_BOX,
                                   composition_mass_fraction={str(z): f for z, f in CONCRETE}),
             energy_grid_MeV=list(D.EGRID), energy_prob=[float(x) for x in d.pE],
             tally_index={'F304': labels, 'F314': labels, 'F324': list(SPECTRUM_DETECTORS)}, **source_meta(d))
    return d, m


def check_p0():
    """rebuild P0 and compare it with the official deck"""
    d = scenes.deck_P0(1e9, 300, False)
    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, 'P0.i')
        d.write(fn)
        with open(fn, 'rb') as f:
            new = f.read()
    with open(paths.INPUTS / 'P0_paciente_em_pe.i', 'rb') as f:
        return new == f.read()


def validate_decks(folder):
    """deck checks of validate.py; the lattice file is read from 1_MCNP_inputs"""
    import validate
    orig = validate.read_cards
    validate.read_cards = lambda fn, base: orig(fn, str(paths.INPUTS))
    for fn in DECKS.values():
        r = validate.check_deck(os.path.join(folder, fn + '.i'))
        print(f'{fn}.i', 'OK' if not r['msgs'] else 'PROBLEM',
              {k: r[k] for k in ('ncells', 'nsurf', 'nent', 'ntally', 'npts')}, r['msgs'])


def main(out):
    same = check_p0()
    print('P0 rebuilt identical to the official deck:', 'yes' if same else 'NO')
    if not same:
        sys.exit('P0 is not reproduced; the test decks were not written.')
    os.makedirs(out, exist_ok=True)
    meta = {}
    for key, fn in DECKS.items():
        d, m = deck_R1() if key == 'R1' else deck_R2()
        d.write(os.path.join(out, fn + '.i'))
        m.update(nps=NPS, ctme_min=CTME, deck=fn + '.i')
        meta[fn] = m
        st = m['source_stats']
        print(f'{fn}.i written: {len(m["regions"])} regions, {m["n_source_entries"]} source entries, '
              f'{st["element_trials"]:.0f} element trials per photon')
    with open(os.path.join(out, META_FILE), 'w', encoding='utf-8', newline='\n') as f:
        json.dump(meta, f, indent=1, ensure_ascii=False, default=float)
    print(f'{META_FILE} written')
    validate_decks(out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else str(paths.DECKS_OUT))
