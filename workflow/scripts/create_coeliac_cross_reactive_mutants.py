#! /usr/bin/env python
"""Create all possible mutants for the cross-reactive epitopes."""

import itertools
import logging

import pandas as pd
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading structure summary from %s', snakemake.input.structures_summary)
structures = pd.read_csv(snakemake.input.structures_summary)
structures.index = structures['path'].map(lambda path: path.rsplit('.', 1)[0])

cdr3_alpha = structures.loc['DQ2_5a1_572G2_DECAB']['cdr3_alpha']
cdr3_beta = structures.loc['DQ2_5a1_572G2_DECAB']['cdr3_beta']

logger.info('Loading sequences from %s', snakemake.input.coeliac_sequences)
sequences = pd.read_csv(snakemake.input.coeliac_sequences, index_col='id')

mhc_chain_1 = sequences['mhc_chain_1'].unique()[0]
mhc_chain_2 = sequences['mhc_chain_2'].unique()[0]

peptides = sequences[['pmhc_calls', 'peptide']].drop_duplicates().set_index('pmhc_calls')['peptide'].to_dict()

tcr_alpha_chain = sequences.loc['coeliac_643_1']['tcr_alpha_chain']
tcr_beta_chain = sequences.loc['coeliac_643_1']['tcr_beta_chain']

logger.debug('Modifying sequences with CDR3a from structure')
start = tcr_alpha_chain.find(sequences.loc['coeliac_643_1']['cdr3a'])
end = start + len(sequences.loc['coeliac_643_1']['cdr3a'])

tcr_alpha_chain = tcr_alpha_chain[: start + 1] + cdr3_alpha + tcr_alpha_chain[end - 1 :]

logger.info('Creating mutations')
differing_positions = [
    (i, (residue_1, residue_2))
    for i, (residue_1, residue_2) in enumerate(zip(peptides['a1'], peptides['w1'], strict=True))
    if residue_1 != residue_2
]

peptide_mutations = []

for peptide in peptides['a1'], peptides['w1']:
    for i, residues in differing_positions:
        for amino_acid in IUPACData.protein_letters:
            if amino_acid in residues:
                continue

            new_peptide = peptide[:i] + amino_acid + peptide[i + 1 :]
            peptide_mutations.append(new_peptide)

peptide_mutations = [
    peptides['a1'][:7] + option[0] + peptides['a1'][8] + option[1] + peptides['a1'][10:]
    for option in itertools.combinations_with_replacement(IUPACData.protein_letters, len(differing_positions))
    if option not in (('L', 'Y'), ('Q', 'F'))
]

logger.info('Exporting mutations to %s', snakemake.output[0])
mutations = pd.DataFrame(
    {
        'id': [f'coeliac_peptide_mut_{num}' for num, _ in enumerate(peptide_mutations, 1)],
        'tcr_alpha_chain': [tcr_alpha_chain for _ in peptide_mutations],
        'tcr_beta_chain': [tcr_beta_chain for _ in peptide_mutations],
        'peptide': peptide_mutations,
        'mhc_chain_1': [mhc_chain_1 for _ in peptide_mutations],
        'mhc_chain_2': [mhc_chain_2 for _ in peptide_mutations],
        'mhc_type': ['MH2' for _ in peptide_mutations],
    }
)
mutations.to_csv(snakemake.output[0], index=False)
