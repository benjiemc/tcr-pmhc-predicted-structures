#!/usr/bin/env python
"""Calculate the TCR:pMHC shape complementarity and buried surface area."""

import contextlib
import json
import logging
import os
import random

import anarci
import numpy as np
from Bio.PDB import Model, PDBParser, Structure
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import IUPACData
from scasa.scasa import Complex
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])


def calculate_bsa(structure: Structure.Structure, chain_annotations: dict[str, str]) -> float:
    """Calculate the buried surface area of the TCR:pMHC interaface."""
    tcr = Structure.Structure('tcr')
    tcr_model = Model.Model(0)

    pmhc = Structure.Structure('pMHC')
    pmhc_model = Model.Model(0)

    for chain in structure[0]:
        if chain_annotations[chain.id] in ('A', 'B', 'G', 'D'):
            tcr_model.add(chain.copy())

        elif chain_annotations[chain.id] in ('MHC', 'antigen'):
            pmhc_model.add(chain.copy())

    tcr.add(tcr_model)
    pmhc.add(pmhc_model)

    sr = ShrakeRupley()

    try:
        sr.compute(structure, level='S')
        sr.compute(tcr, level='S')
        sr.compute(pmhc, level='S')

    except ValueError:
        logger.warning('TCR:pMHC interaction could not be established. Setting buried surface area to 0.0.')
        bsa = 0.0

    else:
        tcr_sasa = tcr.sasa
        overall_sasa = structure.sasa
        pmhc_sasa = pmhc.sasa
        bsa = (tcr_sasa + pmhc_sasa) - overall_sasa

    return bsa


def calculate_sc(path: str, chain_annotations: dict[str, str], seed: int | None = None) -> float:
    """Calculate the shape complementarity of the TCR:pMHC interface using https://github.com/t-whalley/SCASA."""
    if seed:
        logger.debug('Setting random and numpy seeds to %d', seed)
        random.seed(seed)
        np.random.seed(seed)  # noqa: NPY002 - need to use seed because there is no exposed random object

    tcr_pmhc_complex = Complex(
        pdb_file=path,
        complex_1=''.join([chain_id for chain_id, chain_type in chain_annotations.items() if chain_type in ('A', 'B')]),
        complex_2=''.join(
            [chain_id for chain_id, chain_type in chain_annotations.items() if chain_type in ('MHC', 'antigen')]
        ),
    )
    with open(os.devnull, 'w') as fh, contextlib.redirect_stdout(fh):
        try:
            sc = tcr_pmhc_complex.sc()

        except UnboundLocalError:
            logger.warning('TCR:pMHC interaction could not be established. Setting shape complementarity to NaN.')
            sc = float('nan')

    return sc


logger.info('Loading structure from %s', snakemake.input[0])
parser = PDBParser(QUIET=True)
structure = parser.get_structure('', snakemake.input[0])

logger.info('Annotating chains using ANARCI')
chain_annotations = {}
for chain in structure[0]:
    sequence = ''.join([IUPACData.protein_letters_3to1[residue.get_resname().title()] for residue in chain])
    numbering, chain_type = anarci.number(sequence)

    if chain_type in ('MH1', 'CD1', 'B2M', 'GA', 'GB', 'GA1L', 'GA2L', 'GA1', 'GA2', 'MR1'):
        chain_type = 'MHC'

    if not numbering and not chain_type:
        chain_type = 'antigen'

    chain_annotations[chain.id] = chain_type

logger.debug(
    'Chain annotations: %s',
    ' '.join(f'{chain_id}={chain_type}' for chain_id, chain_type in sorted(chain_annotations.items())),
)

logger.info('Calculating shape-complementarity')
sc = calculate_sc(snakemake.input[0], chain_annotations, seed=snakemake.params.seed)

logger.info('Calculating buried surface area')
bsa = calculate_bsa(structure, chain_annotations)

logger.info('Outputing to %s', snakemake.output[0])
with open(snakemake.output[0], 'w') as fh:
    json.dump({'bsa': bsa, 'sc': sc}, fh)
