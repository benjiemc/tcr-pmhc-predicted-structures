#!/usr/bin/env python
"""Add negatives to dash dataset."""

import logging

import anarci
import numpy as np
import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3

logger = logging.getLogger()
setup_logger(logger, level='info')


def collate_cdrs(sequence: str) -> str:
    """Collate CDRs based on TCR alpha chain or beta chain sequence."""
    numbering, _ = anarci.number(sequence)

    cdr1 = []
    cdr2 = []
    cdr3 = []
    for (seq_id, _), amino_acid in numbering:
        if amino_acid == '-':
            continue

        if seq_id in IMGT_CDR1:
            cdr1.append(amino_acid)

        elif seq_id in IMGT_CDR2:
            cdr2.append(amino_acid)

        elif seq_id in IMGT_CDR3:
            cdr3.append(amino_acid)

    return f"{''.join(cdr1)}-{''.join(cdr2)}-{''.join(cdr3)}"


logger.info('Loading Dash et al. sequences from %s', snakemake.input.dash_sequences)
dash_sequences = pd.read_csv(snakemake.input.dash_sequences, index_col='id')

dash_sequences['collated_cdrs'] = dash_sequences.apply(
    lambda row: f'{collate_cdrs(row.tcr_alpha_chain)}-{collate_cdrs(row.tcr_beta_chain)}', axis=1
)
dash_sequences['label'] = 1

dash_pmhc_columns = ['peptide', 'mhc_type', 'mhc_chain_1', 'mhc_chain_2']
dash_peptides = dash_sequences['peptide'].unique()

logger.info('Loading sequences from %s', snakemake.input.sequences)
sequences = pd.read_csv(snakemake.input.sequences, index_col='index')

logger.info('Loading full sequences from %s', snakemake.input.sequences_full)
sequences_full = pd.read_csv(snakemake.input.sequences_full, index_col='id')

logger.info('Selecting positive human full length sequences')
sequences = sequences.merge(
    sequences_full[[column for column in sequences_full.columns if column not in ('peptide', 'mhc_type')]],
    how='inner',
    left_index=True,
    right_index=True,
)

sequences = sequences[sequences['label'].astype(bool)].copy()
sequences = sequences[sequences['species'] == 'Human'].copy()

tcr_columns = [
    'cdr1_alpha',
    'cdr2_alpha',
    'cdr3_alpha',
    'cdr1_beta',
    'cdr2_beta',
    'cdr3_beta',
    'v_alpha',
    'j_alpha',
    'v_beta',
    'j_beta',
    'tcr_alpha_chain',
    'tcr_beta_chain',
    'species',
    'source',
]

sequences['collated_cdrs'] = (
    sequences['cdr1_alpha']
    + '-'
    + sequences['cdr2_alpha']
    + '-'
    + sequences['cdr3_alpha']
    + '-'
    + sequences['cdr1_beta']
    + '-'
    + sequences['cdr2_beta']
    + '-'
    + sequences['cdr3_beta']
    + '-'
)

logger.info('Creating negatives by sampling full length sequences with seed %d', snakemake.params.seed)
rng = np.random.default_rng(snakemake.params.seed)

negative_data = []
for peptide in dash_peptides:
    in_group = dash_sequences[dash_sequences['peptide'] == peptide]
    out_group = sequences[
        (sequences['peptide'] != peptide) & (~sequences['collated_cdrs'].isin(in_group['collated_cdrs']))
    ]

    negatives = (out_group.sample(len(in_group), random_state=rng) if len(out_group) >= len(in_group) else out_group)[
        tcr_columns
    ]
    negatives = negatives.rename({'species': 'tcr_species'}, axis='columns')
    negatives['id'] = in_group.index + '_' + negatives.index.astype(str)

    negatives = pd.concat(
        [
            negatives,
            pd.DataFrame(
                [in_group[dash_pmhc_columns].iloc[0].to_numpy()] * len(negatives),
                columns=in_group[dash_pmhc_columns].iloc[0].index,
                index=negatives.index,
            ),
        ],
        axis=1,
    )

    negatives = negatives.set_index('id')
    negatives['source'] = 'shuffling'
    negatives['label'] = 0

    negative_data.append(negatives)

negative_data = pd.concat(negative_data)

dash_sequences = dash_sequences.drop('collated_cdrs', axis='columns')

logger.info('Combinding Dash et al. sequences and negatives')
dash_sequences_with_negatives = pd.concat([dash_sequences, negative_data[dash_sequences.columns]])

logger.info('Outputting to %s', snakemake.output[0])
dash_sequences_with_negatives.to_csv(snakemake.output[0])
