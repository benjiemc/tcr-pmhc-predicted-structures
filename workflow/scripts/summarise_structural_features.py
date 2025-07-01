#!/usr/bin/env python
"""Calculate structural properties of TCR:pMHC structures."""

import json
import warnings

import stcrpy
from snakemake.script import snakemake

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    tcr = stcrpy.load_TCR(snakemake.input[0])

tcr_info = {}

docking_geometry = tcr.calculate_docking_geometry(mode='cys')

tcr_info['pitch_angle'] = float(docking_geometry['pitch_angle'])
tcr_info['scanning_angle'] = float(docking_geometry['scanning_angle'])
tcr_info['z_score'] = float(docking_geometry['tcr_com'][0][2])

with open(snakemake.output[0], 'w') as fh:
    json.dump(tcr_info, fh)
