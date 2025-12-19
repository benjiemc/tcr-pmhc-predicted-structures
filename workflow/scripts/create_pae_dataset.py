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
sequences = pd.read_csv(snakemake.input.sequences, index_col='index')
sequences.index = sequences.index.astype('str')

logger.info('Loading PAEs')
paes = []
folds = []
labels = []
for path in snakemake.input.confidences:
    logger.debug('Loading PAE from %s', path)
    basename = os.path.basename(path)

    pae = np.loadtxt(path)
    pae = pae.flatten()

    ids = re.search(r'^(?P<name>[A-Za-z0-9_]+)_(?P<seed>\d+)_(?P<sample>\d+)', basename).groupdict()

    logger.debug('Collecting info from %s in sequences', ids['name'])
    info = sequences.loc[ids['name']]

    paes.append(pae)
    folds.append(info['fold'])
    labels.append(info['label'])

paes = np.stack(paes)
folds = np.array(folds)
labels = np.array(labels)

logger.info('Inverting and normalising the PAEs')
max_value = np.nanmax(paes)
paes = (max_value - paes) / max_value

logger.info('Filling NaNs as zeros')
paes = np.nan_to_num(paes)

logger.info('Writing dataset to %s', snakemake.output[0])
with h5py.File(snakemake.output[0], 'w') as fh:
    fh['pae'] = paes
    fh['fold'] = folds
    fh['label'] = labels
