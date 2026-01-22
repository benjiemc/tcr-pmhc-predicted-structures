#!/usr/bin/env python
"""Combine sequences with the AF-Training sequences."""

import logging

import numpy as np
import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.utils import find_common_groups

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading raw sequences from %s', snakemake.input.raw_sequences)
sequences = pd.read_csv(snakemake.input.raw_sequences, index_col='index')
sequences.index = sequences.index.astype(str)
sequences['fold'] = sequences['fold'].astype(str)

logger.info('Creating map between pMHCs and folds')
fold_map = sequences.groupby(['peptide', 'mhc1', 'mhc2']).apply(
    lambda group: group['fold'].unique()[0], include_groups=False
)
fold_map.name = 'fold'
fold_map = fold_map.to_frame().reset_index()

fold_map['mhc1'] = fold_map['mhc1'].str.split(';')
fold_map = fold_map.explode('mhc1')

folds = sorted(fold_map['fold'].unique())

peptide_fold_map = fold_map[['peptide', 'fold']].set_index('peptide')['fold'].to_dict()
mhc1_fold_map = fold_map[['mhc1', 'fold']].set_index('mhc1')['fold'].to_dict()

logger.info('Loading sampled sequences from %s', snakemake.input.sequences_sample)
sequences_sample = pd.read_csv(snakemake.input.sequences_sample, index_col='index')

sequences_sample.index = sequences_sample.index.astype(str)
sequences_sample.index.name = 'id'

sequences_sample['fold'] = sequences_sample['fold'].astype(str)
sequences_sample['source'] = 'Sequences'

logger.info('Loading AF-Training sequences from %s', snakemake.input.af_training_sequences)
sequences_af_training = pd.read_csv(snakemake.input.af_training_sequences, index_col='id')
sequences_af_training['source'] = 'AF Training'

logger.info('Adding cross validation fold information based on raw sequences')
sequences_af_training['mhc_fold'] = sequences_af_training['mhc1'].map(mhc1_fold_map).fillna('')
sequences_af_training['peptide_fold'] = sequences_af_training['peptide'].map(peptide_fold_map).fillna('')

sequences_af_training['fold'] = sequences_af_training[['mhc_fold', 'peptide_fold']].apply(set, axis=1).apply(''.join)
sequences_af_training = sequences_af_training.drop(['mhc_fold', 'peptide_fold'], axis='columns')

fold_counts = {fold: count for fold, count in sequences_af_training.value_counts('fold').to_dict().items() if fold}
for fold in folds:
    if fold not in fold_counts:
        fold_counts[fold] = 0

logger.info('Adding new sequences to cross-validation folds')
new_groups = find_common_groups(sequences_af_training[sequences_af_training['fold'] == ''], ['peptide', 'mhc1'])
group_ids, counts = np.unique(new_groups, return_counts=True)

new_fold_counts = fold_counts.copy()

group_mapping = {}
for group_id, count in zip(group_ids, counts, strict=True):
    min_value = np.inf
    min_fold = None
    for fold, fold_count in new_fold_counts.items():
        if fold_count < min_value:
            min_fold = fold
            min_value = fold_count

    group_mapping[group_id] = min_fold
    new_fold_counts[min_fold] += count

sequences_af_training.loc[sequences_af_training['fold'] == '', 'fold'] = [
    group_mapping[group_id] for group_id in new_groups
]

logger.info('Combining sequence datasets')
sequences_combined = pd.concat([sequences_sample, sequences_af_training])

logger.info('Writing combined dataset to %s', snakemake.output[0])
sequences_combined.to_csv(snakemake.output[0])
