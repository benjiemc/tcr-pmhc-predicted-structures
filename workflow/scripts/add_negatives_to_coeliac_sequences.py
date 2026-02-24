#!/usr/bin/env python
"""Add negatives to coeliac data by shuffling binders from other epitopes."""

import logging

import numpy as np
import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading sequences from %s', snakemake.input[0])
sequences = pd.read_csv(snakemake.input[0], index_col='id')

tcr_columns = [
    'clone_id',
    'donor_id',
    'va_gene',
    'cdr3a',
    'ja_gene',
    'vb_gene',
    'cdr3b',
    'jb_gene',
    'tcr_alpha_chain',
    'tcr_beta_chain',
]
tcr_definition_columns = [
    'va_gene',
    'cdr3a',
    'ja_gene',
    'vb_gene',
    'cdr3b',
    'jb_gene',
]

pmhc_columns = ['pmhc_calls', 'mhc1', 'mhc2', 'mhc_type', 'peptide', 'mhc_chain_1', 'mhc_chain_2']

sequences['tcr_id'] = sequences[tcr_definition_columns].apply('_'.join, axis=1)

logger.info('Creating negatives by shuffling (seed: %d)', snakemake.params.seed)
rng = np.random.default_rng(snakemake.params.seed)

negatives = []
for pmhc_call, group in sequences.groupby('pmhc_calls'):
    out_group_tcrs = sequences.loc[~sequences['tcr_id'].isin(group['tcr_id']), tcr_columns]
    sample_size = min(len(group), len(out_group_tcrs))

    negative_tcrs = out_group_tcrs.sample(sample_size, random_state=rng)
    negative_tcrs.index = negative_tcrs.index + f'_neg_{pmhc_call}'

    negative_pmhcs = pd.DataFrame(
        [group[pmhc_columns].iloc[0].to_numpy()] * sample_size,
        columns=pmhc_columns,
        index=negative_tcrs.index,
    )

    negatives.append(pd.concat([negative_tcrs, negative_pmhcs], axis=1))

negatives = pd.concat(negatives)

sequences = sequences.drop('tcr_id', axis='columns')

sequences['label'] = 1
negatives['label'] = 0


sequences_with_negatives = pd.concat([sequences, negatives])

logger.info('Outputting sequences to %s', snakemake.output[0])
sequences_with_negatives.to_csv(snakemake.output[0])
