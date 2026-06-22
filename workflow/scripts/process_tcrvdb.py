#!/usr/bin/env python
"""Process TCRvdb data to follow project conventions."""

import logging

import pandas as pd
import tidytcells
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.utils import assign_species, get_cdr_sequences, mhc_slug_to_code, stitch_sequence

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

SELECTED_EPITOPES = {'ylq': 'YLQPRTFLL', 'glc': 'GLCTLVAML'}
P_VAL_CUTOFF = 1e-5
FOLD_CHANGE_CUTOFF = 0

logger.info('Loading TCRvdb from %s', snakemake.input.tcrvdb)
tcrvdb = pd.read_csv(snakemake.input.tcrvdb, index_col=0)
tcrvdb = tcrvdb.set_index('name')
tcrvdb.index.name = 'id'

logger.info('Selecting peptides: %s', ' and '.join(SELECTED_EPITOPES.values()))
tcrvdb = tcrvdb[tcrvdb['epitope_aa'].isin(SELECTED_EPITOPES.values())].copy()

logger.info('Standardising gene codes with tidytcells')
tcrvdb['v_gene_alpha'] = tcrvdb['TRAV_IMGT'].map(tidytcells.tr.standardise)
tcrvdb['j_gene_alpha'] = tcrvdb['TRAJ_IMGT'].map(tidytcells.tr.standardise)
tcrvdb['v_gene_beta'] = tcrvdb['TRBV_IMGT'].map(tidytcells.tr.standardise)
tcrvdb['j_gene_beta'] = tcrvdb['TRBJ_IMGT'].map(tidytcells.tr.standardise)

tcrvdb['mhc_gene'] = tcrvdb['hla_long'].map(tidytcells.mh.standardise)

logger.info('Creating full length sequences with Stitchr')
tcrvdb['tcr_alpha_chain'] = tcrvdb.apply(
    lambda row: stitch_sequence(row['v_gene_alpha'], row['j_gene_alpha'], row['cdr3_alpha_aa'], 'human'), axis=1
)
tcrvdb['tcr_beta_chain'] = tcrvdb.apply(
    lambda row: stitch_sequence(row['v_gene_beta'], row['j_gene_beta'], row['cdr3_beta_aa'], 'human'), axis=1
)

logger.info('Extracting CDRs')
tcrvdb[['cdr1_alpha', 'cdr2_alpha', 'cdr3_alpha']] = tcrvdb['tcr_alpha_chain'].map(get_cdr_sequences).apply(pd.Series)
tcrvdb[['cdr1_beta', 'cdr2_beta', 'cdr3_beta']] = tcrvdb['tcr_beta_chain'].map(get_cdr_sequences).apply(pd.Series)

tcrvdb = tcrvdb.rename({'epitope_aa': 'peptide'}, axis='columns')

logger.info('Loading MHC sequences')
mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in snakemake.input.mhc_sequences])

mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Human', 'mhc_code'
].apply(tidytcells.mh.standardise, species='homosapiens')

mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Mouse', 'mhc_code'
].apply(tidytcells.mh.standardise, species='musmusculus')

logger.info('Annotating MHC sequence information')
tcrvdb['mhc_chain_1'] = tcrvdb['mhc_gene'].map(mhc_sequences.set_index('mhc_code')['canonical_sequence'].to_dict())

b2m_sequence = 'IQRTPKIQVYSRHPADIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM'
tcrvdb['mhc_chain_2'] = b2m_sequence

logger.info('Annotating validation status')
tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['ylq'], 'status'] = (
    (tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['ylq'], 'log2FoldChange'] < FOLD_CHANGE_CUTOFF)
    & (tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['ylq'], 'padj'] < P_VAL_CUTOFF)
).map(lambda validated: 'Validated' if validated else 'Non-validated')

tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['glc'], 'status'] = (
    (tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['glc'], 'log2FoldChange'] > FOLD_CHANGE_CUTOFF)
    & (tcrvdb.loc[tcrvdb['peptide'] == SELECTED_EPITOPES['glc'], 'padj'] < P_VAL_CUTOFF)
).map(lambda validated: 'Validated' if validated else 'Non-validated')

logger.info('Outputting processed data to %s', snakemake.output[0])
tcrvdb.to_csv(snakemake.output[0])
