#!/usr/bin/env python
"""Calculate structural properties of TCR:pMHC structures."""

import json
import os
import warnings

import stcrpy
from snakemake.script import snakemake

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    tcr = stcrpy.load_TCR(snakemake.input[0])

tcr_info = {}

tcr_info['name'] = snakemake.wildcards.name
tcr_info['seed'] = snakemake.wildcards.seed
tcr_info['sample'] = snakemake.wildcards.sample

docking_geometry = tcr.calculate_docking_geometry(mode='cys')

tcr_info['pitch_angle'] = float(docking_geometry['pitch_angle'])
tcr_info['scanning_angle'] = float(docking_geometry['scanning_angle'])
tcr_info['z_score'] = float(docking_geometry['tcr_com'][0][2])

tcr.standardise_chain_names()
tcr.crop()

with open(snakemake.output.summary, 'w') as fh:
    json.dump(tcr_info, fh)

output_dir = os.path.dirname(snakemake.output.structure)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

tcr.save(snakemake.output.structure)
