# -*- coding: utf-8 -*-
"""Folder locations, all relative to the top folder of the package (the folder above this one).

    1_MCNP_inputs/      decks V, P0, T and R; lattice files; run metadata
    2_MCNP_results/     MCTAL files
    3_model_scripts/    this folder
    4_analysis/         post-processing
    reference_data/     ICRP 110 phantoms and ICRP 116 coefficients (not distributed, see README.md)
    work/               intermediate files written by the scripts (voxel cache, rebuilt decks); can be deleted

The ICRP data can also be kept elsewhere: set the environment variable ICRP_DATA to the folder that holds
ICRP110/ and ICRP116/.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUTS = ROOT / '1_MCNP_inputs'
RESULTS = ROOT / '2_MCNP_results'
SCRIPTS = ROOT / '3_model_scripts'
ANALYSIS = ROOT / '4_analysis'
RUN_METADATA = INPUTS / 'metadados_runs.json'
PHANTOM_MASSES = INPUTS / 'massas_orgaos_fantomas.json'

# ICRP reference data
REFDATA = Path(os.environ.get('ICRP_DATA') or ROOT / 'reference_data')
ICRP110 = REFDATA / 'ICRP110'                  # AM/AM.dat, AM/AM_organs.dat, AM/AM_media.dat, AF/...
ICRP116 = REFDATA / 'ICRP116'
ICRP116_PHOTON_E = ICRP116 / 'Photons' / 'Photons_Effective dose.txt'
ICRP116_PHOTON_E_KERMA = ICRP116 / 'Photons' / 'Photons_Kerma_Effective dose.txt'
ICRP116_SKELETAL = ICRP116 / 'Skeletal fluence to dose response functions'

if not ICRP116_PHOTON_E.is_file():
    sys.exit(f'ICRP 116 data not found ({ICRP116_PHOTON_E}).\n'
             'They are not distributed with this package: see "ICRP data" in README.md.')

# intermediate files
CACHE = ROOT / 'work' / 'cache'                # .npy voxel arrays of the standing and seated phantoms
DECKS_OUT = ROOT / 'work' / 'decks'            # decks rebuilt by build_all.py (to compare with 1_MCNP_inputs)


def mctal(run):
    """MCTAL file of a run"""
    return RESULTS / f'{run}.m'
