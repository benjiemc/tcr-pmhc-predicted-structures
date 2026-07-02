#!/usr/bin/env python
"""Process data from NF9 specific TCR repertoire described at: https://doi.org/10.1172/jci.insight.202235."""

import logging

import pandas as pd
import tidytcells
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.apps.create_full_length_sequences import B2M_SEQUENCES
from tcr_pmhc_predicted_structures.utils import (
    assign_species,
    mhc_slug_to_code,
    stitch_sequence,
)

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Create data frame with TCR sequences')
tcrs = pd.DataFrame(
    [
        ['P1-15', 'KK-008', 'TRAV12-1', 'TRAJ41', 'CVVNMLRNSGYALNF', 'CASSESGGYEQYF', 'TRBV6-1', 'TRBJ2-7'],
        ['P1-1', 'KK-008', 'TRAV12-1', 'TRAJ12', 'CVVNIIMDSSYKILF', 'CASSEGEGYEQYF', 'TRBV6-4', 'TRBJ2-7'],
        ['2-2', 'GV34', 'TRAV12-1', 'TRAJ33', 'CVVNLFDSNYQLIW', 'CASSEGAGYEQYF', 'TRBV2', 'TRBJ2-7'],
        ['5-3', 'GV34', 'TRAV12-1', 'TRAJ28', 'CVVNMLYSGAGSYQLTF', 'CASSESAGYEAFF', 'TRBV6-1', 'TRBJ1-1'],
        ['19-1', 'VKU48', 'TRAV12-1', 'TRAJ20', 'CVVNCLREDDYIKLSF', 'CASSEGRGYEQYF', 'TRBV6-1', 'TRBJ2-7'],
        ['20-2', 'VKU48', 'TRAV4', 'TRAJ9', 'CLVGDIRHTGGFKTIF', 'CASSEDSPLHF', 'TRBV2', 'TRBJ1-6'],
        ['18-10', 'VKU48', 'TRAV12-1', 'TRAJ28', 'CVVNRLQSGAGSYQLTF', 'CASSEAGGYEQYF', 'TRBV2', 'TRBJ2-7'],
        ['10-11', 'GV33', 'TRAV12-1', 'TRAJ41', 'CVVNLLNSNSGYALNF', 'CASSEAGGYEQYF', 'TRBV6-1', 'TRBJ2-7'],
        ['16-2', 'GV36', 'TRAV12-2', 'TRAJ33', 'CAVNGLKDSNYQLIW', 'CASSEGGGYEQYF', 'TRBV6-1', 'TRBJ2-7'],
        ['10-4', 'GV42', 'TRAV12-1', 'TRAJ12', 'CVVNALEDSSYKILF', 'CASSAGQGYEQYF', 'TRBV2', 'TRBJ2-7'],
        ['11-3', 'VKU7', 'TRAV12-1', 'TRAJ33', 'CVVNALRDSNYQLIW', 'CASSGGQGYEQYF', 'TRBV6-1', 'TRBJ2-7'],
        ['17', 'VKU10', 'TRAV10', 'TRAJ57', 'CVVSHLTQGGSEKLVF', 'CASSEGGRAYEQYF', 'TRBV2', 'TRBJ2-7'],
        ['12-3', 'VKU19', 'TRAV12-3', 'TRAJ44', 'CAFTGTASKLTF', 'CASSPELNEQFF', 'TRBV7-8', 'TRBJ2-1'],
    ],
    columns=['id', 'donor', 'v_gene_alpha', 'j_gene_alpha', 'cdr3_alpha', 'cdr3_beta', 'v_gene_beta', 'j_gene_beta'],
)

logger.info('Create peptide mutants')
peptide = 'NYNYLYRLF'
p5_mutants = [peptide[:4] + residue + peptide[5:] for residue in IUPACData.protein_letters]

logger.info('Annotate MHC information')
hla_allele = 'HLA-A*24:02'

tcrs['mhc_alelle'] = hla_allele
tcrs['mhc_type'] = 'MH1'

logger.info('Load MHC sequences from %s', snakemake.input.mhc_sequences)
mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in snakemake.input.mhc_sequences])

mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Human', 'mhc_code'
].apply(tidytcells.mh.standardize, species='homosapiens')

mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Mouse', 'mhc_code'
].apply(tidytcells.mh.standardize, species='musmusculus')

logger.info('Annotating MHC sequences')
tcrs['mhc_chain_1'] = mhc_sequences[mhc_sequences['mhc_code'] == hla_allele].iloc[0]['canonical_sequence']
tcrs['mhc_chain_2'] = B2M_SEQUENCES['Human']

logger.info('Creating full TCR sequences')
tcrs['tcr_alpha_chain'] = tcrs.apply(
    lambda row: stitch_sequence(row['v_gene_alpha'], row['j_gene_alpha'], row['cdr3_alpha'], species='Human'), axis=1
)
tcrs['tcr_beta_chain'] = tcrs.apply(
    lambda row: stitch_sequence(row['v_gene_beta'], row['j_gene_beta'], row['cdr3_beta'], species='Human'), axis=1
)

logger.info('Annotating peptide and mutant information')
tcrs['peptide'] = peptide

for mutant in p5_mutants:
    mutant_aa = mutant[4]
    tcrs[f'peptide_{mutant_aa}'] = mutant

tcrs = tcrs.melt(
    id_vars=[column for column in tcrs.columns if 'peptide' not in column],
    value_vars=[column for column in tcrs.columns if 'peptide' in column],
    var_name='mutant',
    value_name='peptide_sequence',
).rename({'peptide_sequence': 'peptide'}, axis='columns')

logger.info('Sanitising IDs')
tcrs['mutant'] = tcrs['mutant'].map(lambda mutant: mutant.replace('peptide', '').strip('_'))
tcrs.loc[tcrs['mutant'] != '', 'id'] = (
    tcrs.loc[tcrs['mutant'] != '', 'id'] + '_mut_' + tcrs.loc[tcrs['mutant'] != '', 'mutant']
)
tcrs['id'] = 'nf9_' + tcrs['id']
tcrs['id'] = tcrs['id'].str.replace('-', '_')

logger.info('Exporting data to %s', snakemake.output[0])
tcrs.to_csv(snakemake.output[0], index=False)
