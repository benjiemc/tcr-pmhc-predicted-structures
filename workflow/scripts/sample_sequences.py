#!/usr/bin/env python
"""Sample sequences from main sequence file."""

import logging

import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading sequences from %s', snakemake.input.sequences)
sequences = pd.read_csv(snakemake.input.sequences, index_col='index')

logger.info('Loading bad sequences from %s', snakemake.input.bad_sequences)
with open(snakemake.input.bad_sequences, 'r') as fh:
    bad_sequences = [int(line.strip()) for line in fh if line]

logger.info('Removing bad sequences')
sequences = sequences.loc[~sequences.index.isin(bad_sequences)]

if not snakemake.wildcards.unbalance:
    logger.info('Sampling from folds evenly')
    num_folds = len(sequences['fold'].unique())

    folds_sample_size = int(snakemake.wildcards.sample_size) // num_folds
    group_sample_size = folds_sample_size // 2

    sequences_sample = sequences.groupby(['fold', 'label']).sample(
        group_sample_size, random_state=snakemake.params.seed
    )

else:
    logger.info('Randomly sampling')
    group_sample_size = int(snakemake.wildcards.sample_size) // 2

    sequences_sample = sequences.groupby(['label']).sample(group_sample_size, random_state=snakemake.params.seed)

logger.info('Outputting sampled sequences to %s', snakemake.output[0])
sequences_sample.to_csv(snakemake.output[0])
