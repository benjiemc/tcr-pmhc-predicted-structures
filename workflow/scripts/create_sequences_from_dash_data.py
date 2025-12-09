#!/usr/bin/env python
"""Create sequences from Dash data."""

import logging

import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.apps.create_full_length_sequences import B2M_SEQUENCES
from tcr_pmhc_predicted_structures.utils import mhc_code_to_slug, stitch_sequence

logger = logging.getLogger()
setup_logger(logger, 'debug')

logger.info('Loading data sequences from %s', snakemake.input.sequences)
dash = pd.read_csv(snakemake.input.sequences)
dash.index.name = 'id'
dash = dash.reset_index()
dash['id'] = 'dash_' + dash['id'].astype(str)

logger.debug('Dropping Rows with NaNs')
dash = dash.dropna()

logger.info('Stitching full TCR sequences from V gene, J gene, and CDR3')
logger.debug('Stitching full alpha-chain sequences')
dash['tcr_alpha_chain_seq'] = dash.apply(
    lambda row: stitch_sequence(row['v_a_gene'], row['j_a_gene'], row['cdr3_a_aa'], 'human'), axis=1
)

logger.debug('Stitching full beta-chain sequences')
dash['tcr_beta_chain_seq'] = dash.apply(
    lambda row: stitch_sequence(row['v_b_gene'], row['j_b_gene'], row['cdr3_b_aa'], 'human'), axis=1
)

logger.info('Adding epitope sequences')
# Found here: https://www.nature.com/articles/nature22383/tables/1
epitope_info = pd.DataFrame(
    columns=['name', 'species', 'mhc', 'peptide_seq', 'virus'],
    data=[
        ['BMLF', 'human', 'HLA-A*02:01', 'GLCTLVAML', 'EBV'],
        ['M1', 'human', 'HLA-A*02:01', 'GILGFVFTL', 'IAV'],
        ['pp65', 'human', 'HLA-A*02:01', 'NLVPMVATV', 'CMV'],
    ],
)

logger.debug('Loading MHC sequences')
mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in snakemake.input.mhc_sequences])

epitope_info['mhc_slug'] = epitope_info['mhc'].map(mhc_code_to_slug)
epitope_info['mhc_type'] = 'MH1'

epitope_info = epitope_info.merge(
    mhc_sequences[['canonical_sequence']], how='left', left_on='mhc_slug', right_index=True
).rename({'canonical_sequence': 'mhc_chain_1_seq'}, axis='columns')

dash = dash.merge(epitope_info, how='left', left_on='epitope', right_on='name')

dash['mhc_chain_2_seq'] = dash['species'].str.title().map(B2M_SEQUENCES)

logger.info('Outputting sequences to %s', snakemake.output[0])
output_sequences = dash.filter(regex='_seq$|^mhc_type$|^id$')
output_sequences.columns = [column_name.replace('_seq', '') for column_name in output_sequences.columns]

output_sequences.to_csv(snakemake.output[0], index=False)
