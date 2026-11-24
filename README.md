# Monte Carlo dosimetry of therapeutic radioligands

Input decks, tally results and scripts of the study

**Monte Carlo assessment of ambient dose equivalent and effective dose around patients treated with 177Lu, 161Tb, 212Pb and 225Ac radioligands**.

The simulations were run with MCNP 6.3.0 and the ICRP 110 adult voxel phantoms. The patient phantom carries a photon source that is resolved in energy (27 energies) and tagged by source region (24 regions), so one set of runs gives H\*(10), air kerma and effective dose for any of the four radionuclides and for any activity distribution, by folding. With the files below it is possible to follow every number of the article and of its Supplementary Material from the raw tallies to the printed value, to rebuild the input decks and to rerun the simulations.

Section, figure and table numbers in this file refer to the Medical Physics manuscript and its Supplementary Material. The scripts and reference outputs retain their original table labels. Former manuscript Tables 1–4 are now Supplementary Tables S26–S29, and former Table 5 is now manuscript Table 1.

## Contents

| Folder or file | What it holds |
| --- | --- |
| `1_MCNP_inputs/` | The 19 input decks (V01–V08, P0, T01–T08, R1, R2), the run metadata and the batch files used to run them |
| `2_MCNP_results/` | MCTAL file of every run and `run_summary.csv` |
| `3_generator/` | Scripts that define the model (phantoms, source, tallies, scenes) and write the decks, and the scripts that check the decks |
| `4_analysis/` | Post-processing, from the MCTAL files to the tables and figures, with reference outputs |
| `reference_data/` | Place for the ICRP data, which are not distributed here |
| `Supplementary_material.docx` | Supplementary Material of the article (Sections S1–S9, Figs. S1–S3, Tables S1–S29) |

## Requirements

* Python 3.9 or later with NumPy 2.0 or later (the code uses `numpy.trapezoid`); no other package is needed.
* The ICRP data described in the next section.
* MCNP 6.3.0 with the photoatomic library `eprdata14`, only to rerun the simulations (distributed by RSICC and the NEA Data Bank under licence).

## ICRP data

The scripts read two data sets that ICRP publishes as supplementary data and that we do not redistribute:

* Publication 110, adult reference computational phantoms (version 1.2): folders `AM` and `AF`;
* Publication 116, conversion coefficients (version 2): folders `Photons` and `Skeletal fluence to dose response functions`.

Download them from the ICRP website and unpack them into `reference_data/` so that these files exist:

```
reference_data/ICRP110/AM/AM.dat            (and AM_organs.dat, AM_media.dat; the same for AF)
reference_data/ICRP116/Photons/Photons_Effective dose.txt
reference_data/ICRP116/Photons/Photons_Kerma_Effective dose.txt
reference_data/ICRP116/Skeletal fluence to dose response functions/Table D.1 - Regional DRF.txt
reference_data/ICRP116/Skeletal fluence to dose response functions/Table 3.2.txt
```

If the data are kept elsewhere, set the environment variable `ICRP_DATA` to the folder that contains `ICRP110` and `ICRP116`. The post-processing needs only the Publication 116 files (the air-kerma coefficients and the skeletal response functions are read from them). The Publication 110 files are needed to rebuild the decks and the lattice files.

The decks read the voxel geometry from three lattice files (`lat_PAMstd.txt`, `lat_TAMstd.txt` and `lat_TAFstd.txt`, 8 MB in total) with `read file=`. These files are the ICRP 110 voxel arrays rearranged as MCNP lattice cards, so they are not in the repository either. `build_all.py` writes them in about two minutes (see "Rebuilding the input decks").

## Post-processing

From `4_analysis`, in this order:

```
python a1_01_verification_icrp116.py        > my_a1_01.txt
python a1_02_benchmark_and_source_models.py > my_a1_02.txt
python a1_03_article1_results.py            > my_a1_03.txt
python a1_04_supporting_numbers.py          > my_a1_04.txt
python a1_05_supplementary_tests.py         > my_a1_05.txt
python a1_06_supplementary_tables.py        > my_a1_06.txt
```

The whole sequence takes about 15 s on a workstation. `a1_04` reads `article1_results.json`, written by `a1_03`, and `a1_06` reads the files written by `a1_02`–`a1_05`, so the order matters. The output of `a1_06` is UTF-8.

Each output can be compared line by line with the corresponding `output_a1_0x_….txt`. We obtained the same outputs with Python 3.13.17 / NumPy 2.5.3 and Python 3.10.12 / NumPy 2.2.6 on Linux and with Python 3.11.17 / NumPy 2.4.6 on Windows. The only line that depends on the version is the consistency check at the end of `a1_06`, which prints 0 or 2.2e-16. The CSV files are reproduced byte for byte. The values in the JSON files can differ in the last digit (relative difference below 10⁻¹⁵), because the order of the floating-point sums changes with the NumPy version.

The scripts rewrite the JSON and CSV files of `4_analysis`; work on a copy to keep the reference files.

## Rebuilding the input decks

From `3_model_scripts`:

```
python build_all.py
python write_meta.py
python build_tests.py
```

The output goes to `work/decks` (or to a folder given as argument). The first run reads the ICRP 110 phantoms and caches the voxel arrays in `work/cache` (about 0.2 GB). `build_all.py` writes the production decks, the lattice files and `metadados_runs.json`; `write_meta.py` writes `massas_orgaos_fantomas.json`; `build_tests.py` writes the decks of the two supplementary runs and `metadados_testes.json`, after checking that P0 is rebuilt identical to the deck in `1_MCNP_inputs`. The rebuilt decks and lattice files are identical, byte for byte, to those that were run; the JSON files have the same content, and their floating-point values can differ in the last digit with the NumPy version.

To run MCNP or the deck checks, copy the three lattice files to the deck folder:

```
copy ..\work\decks\lat_*std.txt ..\1_MCNP_inputs
```

`build_tests.py` runs the checks of `validate.py` on the two supplementary decks and therefore needs the lattice files in `1_MCNP_inputs`; run it after this copy. The deck checks are run with `python validate.py` (all decks of `1_MCNP_inputs`) or `python check2.py <deck> ...`, and likewise for `check3.py`, `overlap2.py` and `containment.py`. Their results for the decks that were run are in the `validation_*.log` files.

The same scripts also hold the scenes of a companion study on passengers travelling with a patient (runs B01–B07, seated phantoms). `build_all.py` writes those decks too and the logs and `metadados_runs.json` list them; their results are not part of this repository.

## Rerunning the simulations

Run `02_FASE1_verificacao_ICRP116.bat`, `03_FASE2_artigo1_paciente_em_pe.bat` and `07_TESTES_R1_R2.bat` from `1_MCNP_inputs`, after setting the MCNP executable and the number of threads at the top of each file. A run is skipped when its MCTAL file already exists in `2_MCNP_results`, so the results in the repository are never overwritten; move them away to rerun. The full output and the restart files go to `mcnp_output/`.

On 20 threads the 17 production runs took 3 h 16 min and each supplementary run about 1 h. All runs stop at a computer-time limit (CTME), so a rerun on another machine processes a different number of histories and reproduces the results within their statistical uncertainties, not bit for bit.

## Files

### `1_MCNP_inputs`

| File | Description | Article |
| --- | --- | --- |
| `V01_AM_AP.i` … `V08_AF_RLAT.i` | Verification: broad parallel photon beams (AP, PA, LLAT, RLAT) on the AM and AF phantoms, 18 energies from 20 keV to 3 MeV | Sections 2.7 and 3.1; Fig. 3a |
| `P0_paciente_em_pe.i` | Standing AM patient alone in air, with the energy-resolved, region-tagged source and 76 air-sphere detectors | Sections 2.3, 2.6, 3.1–3.3; Table S28; Figs. 3b–c and 4; Figs. S2b and S3 |
| `T01_AM_g030.i` … `T08_AF_g200.i` | Patient and three exposed persons, AM (T01–T04) or AF (T05–T08), in front, beside and behind at 30, 50, 100 and 200 cm | Sections 2.6, 3.3–3.5; Table S29 and manuscript Table 1; Figs. 5–6 |
| `R1_P0_Tb161_direto.i` | Geometry of P0; half of the histories emit the 21 photon lines of ¹⁶¹Tb and the other half the same lines folded on the 27-energy grid | Sections 2.7 and 3.1; Section S3.3 |
| `R2_P0_piso_concreto.i` | P0 with a 30 cm concrete floor (8.4 m × 8.4 m) under the patient; a second tally bin separates the photons that have left the floor | Sections 2.7, 3.1 and 4; Section S3.4 |
| `metadados_runs.json` | Tally index (detector labels, organ targets), source regions and probabilities, energy grid and sampling probabilities, region time tags, NPS and CTME, phantom placements | needed to read the MCTAL bins |
| `metadados_testes.json` | The same for R1 and R2, with the floor geometry and composition | needed to read the MCTAL bins of R1 and R2 |
| `massas_orgaos_fantomas.json` | Organ and tissue masses and slab structure of the phantoms, skeletal masses of ICRP 116 | needed for the organ doses (Section 2.5) |
| `*.bat` | Windows batch files that run the decks | – |

### `2_MCNP_results`

One MCTAL file per run (`.m`), with every tally bin and its relative error, and `run_summary.csv`: histories, computer time, termination, start and end, warnings, lost particles, outcome of the ten statistical checks of MCNP for each tally, and the SHA-256 sum of the MCTAL file.

All 19 runs ended at the computer-time limit and lost no particles. Most tallies that missed one of the ten checks missed the trend check of the figure of merit, which contains the computer time and so drifts when the computing speed changes during a run; `run_summary.csv` lists every case. R1 processed 8.6 × 10⁸ histories and MCNP warned that the period of the random-number generator was exceeded (default stride). Each history used 424 random numbers on average, far fewer than the stride, so the later histories start at positions of the sequence that the earlier ones did not reach; the direct-to-folded ratios of R1 agree with unity within their uncertainties (Section 3.1).

### `3_model_scripts`

| File | Description |
| --- | --- |
| `mcnpgen.py` | Deck writer: phantoms as stacks of cropped z-slab lattices, materials, source (SDEF with the energy grid and the region time tags), tallies |
| `scenes.py` | The scenes: verification beams, patient alone, patient with exposed persons |
| `phantoms.py`, `slabs.py`, `phantom_cache.py` | Reading of the ICRP 110 phantoms, partition in slabs, cache of the voxel arrays |
| `data.py` | Photon emissions (DDEP; ¹⁶¹Tb from Westerbergh et al., 2025), fluence-to-H\*(10) coefficients of ICRP 74, air-kerma coefficients and skeletal response functions of ICRP 116, source regions, energy grid and sampling probabilities |
| `build_all.py`, `write_meta.py`, `build_tests.py` | Write the decks, the lattice files and the JSON files |
| `build_seated.py`, `joints.py`, `check_seated.py`, `AM_seated.json`, `AF_seated.json` | Seated phantoms of the companion study |
| `validate.py`, `check2.py`, `check3.py`, `overlap.py`, `overlap2.py`, `containment.py` | Checks of the decks: structure, source definition (sampling boxes, uniformity per organ, region weights), energy grid, region tags and tally sizes, overlap of all objects, containment of the slabs |
| `validation_*.log` | Results of these checks for the decks that were run: no errors and no overlaps |
| `paths.py` | Folder locations |

### `4_analysis`

| File | Description | Article |
| --- | --- | --- |
| `mctal.py` | MCTAL reader | – |
| `postproc_core.py` | Response matrices R(region, energy) per emitted photon, folding with the emission lines, organ doses and sex-averaged effective dose (ICRP 103), reference [¹⁷⁷Lu]Lu-PSMA-617 biokinetic model | Sections 2.3–2.5 |
| `a1_01_verification_icrp116.py` | E/Φ of the V runs relative to ICRP 116 | Section 3.1; Fig. 3a |
| `a1_02_benchmark_and_source_models.py` | H\*(10) of the model patient versus time; H\*(10) around the patient; in-body point and line sources and bare point source; writes `article1_detector_results.json` and `article1_photon_spectra_1m.csv` | Sections 2.4, 3.1–3.3; Tables S27–S28; Figs. 3b–c and 4; Figs. S2b and S3 |
| `a1_03_article1_results.py` | Effective dose of the exposed persons, air kerma and H\*(10)/Ka, point-source constants, H\*(10) versus time, time-integrated doses, release criterion; writes `article1_results.json` | Sections 3.1, 3.2, 3.4, 3.5; Tables S26 and S29 and manuscript Table 1; Figs. 3b–c, 5, 6 |
| `a1_04_supporting_numbers.py` | Numbers quoted in the text: height and direction dependence, simplified sources, H\*(10)/E, cumulative dose fraction; writes `article1_cumulative_fraction.csv` | Sections 3.2–3.5; Figs. 5a, 6b |
| `a1_05_supplementary_tests.py` | Runs R1 (direct versus folded ¹⁶¹Tb spectrum) and R2 (concrete floor); writes `article1_supplementary_tests.csv` | Sections 2.7, 3.1 and 4; Sections S3.3 and S3.4 |
| `a1_06_supplementary_tables.py` | Prints Tables S1–S25 of the Supplementary Material in Markdown and checks the recomputed effective doses against `article1_results.json` | Tables S1–S25 |
| `article1_*.json`, `article1_*.csv` | Results written by the scripts above | – |
| `output_a1_0x_….txt` | Reference output of each script | – |

Keys of `article1_results.json` (rates per GBq in the patient's body with the activity distribution at 4 h and the bladder empty, in µSv h⁻¹ or µGy h⁻¹; integrals in µSv for 7.4 GBq of ¹⁷⁷Lu or ¹⁶¹Tb, 0.2 GBq of ²¹²Pb and 8 MBq of ²²⁵Ac):

* `E["<POSITION>_<gap cm>_<nuclide>"]`: effective dose with the spongiosa method (`spong`), its relative error (`relerr`) and with the skeletal response functions (`drf`); positions FRONT, SIDE and BACK (Table S29);
* `Ka["<detector>_<nuclide>"]`: `Ka`, `H` (H\*(10)) and their `ratio`;
* `Gamma[<nuclide>]`: H\*(10) and air-kerma rate constants of a bare point source at 1 m (Table S26);
* `timecurve`: H\*(10) for 7.4 GBq of ¹⁷⁷Lu versus time at 0.5, 1 and 2 m, bladder empty and full (Fig. 3b–c);
* `integrals["<nuclide>_<release h>"]`: time-integrated E at 30 cm, 50 cm and 1 m, H\*(10) at 1 m and the point-source estimate `H_point_phys` (manuscript Table 1; Fig. 6a);
* `fractions[<nuclide>]`: share of the 30-day E delivered in the first day and week;
* `release2m[<nuclide>]`: H\*(10) at 2 m per GBq in the most exposed direction and the activity that gives 30 µSv h⁻¹ (Section 3.5).

Detector labels read `<front|back|left|right>_g<gap cm>_z<height cm>` or `diag<angle>_R<distance from the trunk axis>_z<height>`; the coordinates are in `metadados_runs.json` and in `article1_detector_results.json`.

## Data sources

| Data | Source | Here? |
| --- | --- | --- |
| Adult reference voxel phantoms | ICRP Publication 110, supplementary data v1.2 | No (see "ICRP data") |
| Effective dose and air-kerma coefficients, skeletal response functions | ICRP Publication 116, supplementary data v2 | No (see "ICRP data") |
| Fluence-to-H\*(10) coefficients | ICRP Publication 74, Table A.21 | In `data.py` |
| Photon emissions | DDEP evaluations (LNHB); ¹⁶¹Tb: Westerbergh et al. (2025) | In `data.py` |
| Biokinetic model of [¹⁷⁷Lu]Lu-PSMA-617 | Built from seven published dosimetry studies (Table S27; Section S2.4) | In `postproc_core.py` (`lu_model`) |
| Measured dose rates used as benchmark | Published values, cited in the article | Quoted in the scripts |

No patient data were used.

## Not included

* The full MCNP output files (`.o`, 0.1 GB) and the restart files (2.1 GB). `run_summary.csv` gives what is usually needed from them; they are available on request.
* The scripts that draw the figures. They only plot the files of `4_analysis`.

## Licence

The scripts are released under the MIT licence (see `LICENSE`). The input decks, the results and the Supplementary Material are released under CC BY 4.0. The ICRP data keep their own terms of use.

## Citation

Please cite the article. The full reference will be added here when it is published.
