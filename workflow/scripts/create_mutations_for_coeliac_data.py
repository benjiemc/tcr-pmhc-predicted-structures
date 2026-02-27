#!/usr/bin/env python
"""Create mutations for coeliac data based off of similar analysis done in the paper."""

import logging

import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading structures summary from %s', snakemake.input.structures_summary)
structures = pd.read_csv(snakemake.input.structures_summary)
structures.index = structures['path'].map(lambda path: path.rsplit('.', 1)[0])

cdr3_alpha = structures.loc['DQ2_5a1_572G2_DECAB']['cdr3_alpha']
cdr3_beta = structures.loc['DQ2_5a1_572G2_DECAB']['cdr3_beta']

logger.info('Loading sequences from %s', snakemake.input.coeliac_sequences)
sequences = pd.read_csv(snakemake.input.coeliac_sequences, index_col='id')

mhc_chain_1 = sequences['mhc_chain_1'].unique()[0]
mhc_chain_2 = sequences['mhc_chain_2'].unique()[0]

peptides = sequences[['pmhc_calls', 'peptide']].drop_duplicates().set_index('pmhc_calls')['peptide'].to_dict()
selected_peptides = {name: peptide for name, peptide in peptides.items() if name in ('a1', 'w1')}

tcr_alpha_chain = sequences.loc['coeliac_643_1']['tcr_alpha_chain']
tcr_beta_chain = sequences.loc['coeliac_643_1']['tcr_beta_chain']

logger.info('Getting new alpha chain sequence by modifying CDR3a')
start = tcr_alpha_chain.find(sequences.loc['coeliac_643_1']['cdr3a'])
end = start + len(sequences.loc['coeliac_643_1']['cdr3a'])

tcr_alpha_chain = tcr_alpha_chain[: start + 1] + cdr3_alpha + tcr_alpha_chain[end - 1 :]

logger.info('Creating alanine mutations')
mutations = {
    'alpha': {'cdr3': [], 'sequence': []},
    'beta': {'cdr3': [], 'sequence': []},
}
for chain_type, starting_cdr3, sequence in (
    ('alpha', cdr3_alpha, tcr_alpha_chain),
    ('beta', cdr3_beta, tcr_beta_chain),
):
    cdr3_start = sequence.find(starting_cdr3)
    cdr3_end = cdr3_start + len(starting_cdr3)

    for i, amino_acid in enumerate(starting_cdr3[1:-1], 1):
        if amino_acid == 'A':
            continue

        new_cdr3 = starting_cdr3[:i] + 'A' + starting_cdr3[i + 1 :]
        new_sequence = sequence[:cdr3_start] + new_cdr3 + sequence[cdr3_end:]

        mutations[chain_type]['cdr3'].append(new_cdr3)
        mutations[chain_type]['sequence'].append(new_sequence)

logger.info('Compiling mutations into ternary complexes')
df = {
    'id': [],
    'pmhc_calls': [],
    'cdr3_alpha': [],
    'cdr3_beta': [],
    'tcr_alpha_chain': [],
    'tcr_beta_chain': [],
    'peptide': [],
    'mhc_chain_1': [],
    'mhc_chain_2': [],
    'mhc_type': [],
}
for pmhc_calls, peptide in selected_peptides.items():
    df['id'].append(f'coeliac_{pmhc_calls}_mut_original')

    df['cdr3_alpha'].append(cdr3_alpha)
    df['cdr3_beta'].append(cdr3_beta)
    df['pmhc_calls'].append(pmhc_calls)

    df['tcr_alpha_chain'].append(tcr_alpha_chain)
    df['tcr_beta_chain'].append(tcr_beta_chain)
    df['peptide'].append(peptide)
    df['mhc_chain_1'].append(mhc_chain_1)
    df['mhc_chain_2'].append(mhc_chain_2)
    df['mhc_type'].append('MH2')

    count = 1
    for chain_type, entries in mutations.items():
        for cdr3, sequence in zip(entries['cdr3'], entries['sequence'], strict=True):
            df['id'].append(f'coeliac_{pmhc_calls}_mut_{count}')
            df['pmhc_calls'].append(pmhc_calls)
            df['peptide'].append(peptide)
            df['mhc_chain_1'].append(mhc_chain_1)
            df['mhc_chain_2'].append(mhc_chain_2)
            df['mhc_type'].append('MH2')

            match chain_type:
                case 'alpha':
                    df['cdr3_alpha'].append(cdr3)
                    df['tcr_alpha_chain'].append(sequence)

                    df['cdr3_beta'].append(cdr3_beta)
                    df['tcr_beta_chain'].append(tcr_beta_chain)

                case 'beta':
                    df['cdr3_alpha'].append(cdr3_alpha)
                    df['tcr_alpha_chain'].append(tcr_alpha_chain)

                    df['cdr3_beta'].append(cdr3)
                    df['tcr_beta_chain'].append(sequence)

            count += 1

df = pd.DataFrame(df)
df = df.set_index('id')

logger.info('Outputting dataset to %s', snakemake.output[0])
df.to_csv(snakemake.output[0])
