#!/usr/bin/env python
"""Create an h5 dataset from PAE matrices."""

import logging
import os
import re

import h5py
import numpy as np
import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level='info')

logger.info('Loading sequences from %s', snakemake.input.sequences)
sequences = pd.read_csv(snakemake.input.sequences, index_col='id')
sequences.index = sequences.index.astype('str')

logger.info('Loading PAEs')
ids = []
seeds = []
samples = []
sources = []
paes = []
folds = []
labels = []
for path in snakemake.input.confidences:
    logger.debug('Loading PAE from %s', path)
    basename = os.path.basename(path)

    pae = np.loadtxt(path)
    pae = pae.flatten()

    entry_name = re.search(r'^(?P<name>[A-Za-z0-9_]+)_(?P<seed>\d+)_(?P<sample>\d+)', basename).groupdict()

    logger.debug('Collecting info from %s in sequences', entry_name['name'])
    info = sequences.loc[entry_name['name']]

    ids.append(entry_name['name'])
    seeds.append(int(entry_name['seed']))
    samples.append(int(entry_name['sample']))

    sources.append(info['source'])
    folds.append(info['fold'])
    labels.append(info['label'])

    paes.append(pae)

seeds = np.array(seeds)
samples = np.array(samples)
folds = np.array(folds)
labels = np.array(labels)

paes = np.stack(paes)

logger.info('Inverting and normalising the PAEs')
max_value = np.nanmax(paes)
paes = (max_value - paes) / max_value

logger.info('Filling NaNs as zeros')
paes = np.nan_to_num(paes)

logger.info('Writing dataset to %s', snakemake.output[0])
with h5py.File(snakemake.output[0], 'w') as fh:
    fh['id'] = ids
    fh['seed'] = seeds
    fh['sample'] = samples
    fh['pae'] = paes
    fh['source'] = sources
    fh['fold'] = folds
    fh['label'] = labels
