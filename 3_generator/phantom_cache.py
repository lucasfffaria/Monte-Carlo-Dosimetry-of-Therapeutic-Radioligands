# -*- coding: utf-8 -*-
"""Voxel arrays of the ICRP 110 phantoms (standing and seated), cached in work/cache.

- standing: read from the original ICRP 110 files (AM.dat, AF.dat), see paths.py;
- seated: built by build_seated.build() with the joint angles recorded in AM_seated.json and AF_seated.json.

The cached .npy files are regenerable intermediates and can be deleted at any time.
"""
import json
import numpy as np
import paths

DIMS = {'AM': (254, 127, 222), 'AF': (299, 137, 348)}      # columns, rows, slices (ICRP 110)


def std_array(sex):
    """array [slice, row, column] of organ identifiers of the standing reference phantom"""
    fn = paths.CACHE / f'{sex}.npy'
    if not fn.exists():
        nc, nr, ns = DIMS[sex]
        a = np.fromfile(paths.ICRP110 / sex / f'{sex}.dat', sep=' ', dtype=np.int32)
        assert a.size == nc * nr * ns, (sex, a.size)
        paths.CACHE.mkdir(parents=True, exist_ok=True)
        np.save(fn, a.reshape(ns, nr, nc))
    return np.load(fn)


def seated_info(sex):
    with open(paths.SCRIPTS / f'{sex}_seated.json', encoding='utf-8') as f:
        return json.load(f)


def seated_array(sex):
    """array of the seated phantom (same convention), rebuilt from the parameters in the .json file"""
    fn = paths.CACHE / f'{sex}_seated.npy'
    if not fn.exists():
        from build_seated import build
        info = seated_info(sex)
        out, _src, new = build(sex, info['th_hip'], info['th_elb'], info['psi'], verbose=False)
        assert list(new['origin_idx']) == list(info['origin_idx']), (sex, new['origin_idx'], info['origin_idx'])
        paths.CACHE.mkdir(parents=True, exist_ok=True)
        np.save(fn, out)
    return np.load(fn)
