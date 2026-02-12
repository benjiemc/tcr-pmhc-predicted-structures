#!/usr/bin/env python
"""Create an h5 dataset from PAE matrices, use all of the PAE matrix and left pad to the max size."""

import json
import logging
import os
import re

import h5py
import numpy as np
import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Getting max size and max value')
num_confidences = len(snakemake.input.confidences)

max_size = 0
max_value = -np.inf
for num, path in enumerate(snakemake.input.confidences, 1):
    logger.debug('Loading confidence scores from %s - %d of %d', path, num, num_confidences)
    with open(path, 'r') as fh:
        confidences = json.load(fh)

    pae = np.array(confidences['pae'])

    max_size = max(max_size, pae.shape[0])
    max_value = max(max_value, np.nanmax(pae))

logger.debug('Max size is %d, max value is %f', max_size, max_value)

logger.info('Loading sequences from %s', snakemake.input.sequences)
sequences = pd.read_csv(snakemake.input.sequences, index_col='id')
sequences.index = sequences.index.astype('str')

logger.info('Creating dataset at %s', snakemake.output[0])
with h5py.File(snakemake.output[0], 'w') as fh:
    str_dt = h5py.string_dtype(encoding='utf-8')
    ids = fh.create_dataset('id', shape=(num_confidences,), dtype=str_dt)
    seeds = fh.create_dataset('seed', shape=(num_confidences,), dtype=np.int32)
    samples = fh.create_dataset('sample', shape=(num_confidences,), dtype=np.int32)
    sources = fh.create_dataset('source', shape=(num_confidences,), dtype=np.int8)
    folds = fh.create_dataset('fold', shape=(num_confidences,), dtype=np.int32)
    labels = fh.create_dataset('label', shape=(num_confidences,), dtype=np.int8)
    paes = fh.create_dataset(
        'pae',
        shape=(num_confidences, max_size * max_size),
        dtype=np.float32,
        chunks=(1, max_size * max_size),
        compression='gzip',
        compression_opts=4,
        fillvalue=0.0,
    )

    for i, path in enumerate(snakemake.input.confidences):
        logger.debug('Adding %d of %d', i + 1, num_confidences)

        entry_name = re.search(
            r'^(?P<name>[A-Za-z0-9_]+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_confidences.json',
            os.path.basename(path),
        ).groupdict()

        logger.debug('Collecting info from %s in sequences', entry_name['name'])
        info = sequences.loc[entry_name['name']]

        logger.debug('Loading confidence scores from %s', path)
        with open(path, 'r') as fh_in:
            confidences = json.load(fh_in)

        pae = np.array(confidences['pae'], dtype=np.float32)

        logger.debug('Padding PAE')
        pae = np.pad(
            pae, ((0, max_size - pae.shape[0]), (0, max_size - pae.shape[0])), mode='constant', constant_values=np.nan
        )

        logger.debug('Inverting and normalising the PAE')
        pae = (max_value - pae) / max_value

        logger.debug('Filling NaNs as zeros')
        pae = np.nan_to_num(pae, nan=0.0)

        logger.debug('Flattening PAE matrices')
        pae = pae.flatten()

        logger.debug('Writing line')
        ids[i] = entry_name['name']
        seeds[i] = int(entry_name['seed'])
        samples[i] = int(entry_name['sample'])

        sources[i] = int(info['source'] == 'Sequences')
        folds[i] = info['fold']
        labels[i] = info['label']

        paes[i, :] = pae
