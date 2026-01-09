#!/usr/bin/env python
"""Calculate the DockQ score between a predicted structure and a reference one."""

import json

import stcrpy
from snakemake.script import snakemake
from stcrpy.tcr_metrics.tcr_dockq import TCRDockQ

tcr = stcrpy.load_TCR(snakemake.input.prediction, include_symmetry_mates=False)
tcr_ref = stcrpy.load_TCR(snakemake.input.structure, include_symmetry_mates=False)

if tcr:
    if tcr.get_antigen() and tcr.get_MHC():
        dockq_result = TCRDockQ().tcr_dockq(tcr, tcr_ref)
        dockq_result = dockq_result['best_result']['TM']

    else:
        print('No antigen or MHC detected')
        dockq_result = {}

else:
    print('TCR not recognised')
    dockq_result = {}

with open(snakemake.output[0], 'w') as fh:
    json.dump(dockq_result, fh)
