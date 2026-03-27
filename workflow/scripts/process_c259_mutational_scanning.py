#!/usr/bin/env python
"""Process mutational scanning for c259 TCR."""

import logging
import re

import pandas as pd
import tidytcells
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.utils import assign_species, get_cdr_sequences, mhc_slug_to_code

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading mutational affinity data from %s', snakemake.input.mutational_scans)
with open(snakemake.input.mutational_scans, 'r') as fh:
    data = [line.strip('\n') for line in fh if line]

logger.info('Processing affinity data')
data_formated = {'position': []}
for line in data:
    if line.startswith('Position'):
        positions = [position.replace(' ', '_').lower() for position in re.split(r'(?<=Position \d)\s', line)]

    elif line.startswith('Abbr'):
        header = line.replace('Abbr', 'id').replace('Geo. ', 'geometric_').lower().split()

        for name in header:
            if name not in data_formated:
                data_formated[name] = []

    else:
        values = line.replace('N/A', 'NaN NaN NaN').split()

        num_entries = len(values)
        num_positions = len(positions)
        period = num_entries // num_positions

        for i, (name, value) in enumerate(zip(header, values, strict=True)):
            try:
                value_converted = float(value)

            except ValueError:
                value_converted = value

            data_formated[name].append(value_converted)

            if i % period == 0:
                data_formated['position'].append(positions[i // period])

data_formated = pd.DataFrame(data_formated)
data_formated = data_formated.sort_values('id')
data_formated = data_formated.set_index('id')

data_formated = data_formated.rename({'n': 'num_samples', 'sequence': 'peptide'}, axis='columns')
data_formated['num_samples'] = data_formated['num_samples'].fillna(0).astype(int)

original = pd.Series(
    ['original_9V', 'SLLMWITQV', 1.240, 0.129, 3],
    index=['id', 'peptide', 'geometric_mean', 'geometric_sd', 'num_samples'],
)

data_formated = pd.concat([original.to_frame().T.set_index('id'), data_formated])

logger.info('Adding sequence information')
data_formated['mhc_type'] = 'MH1'
data_formated['mhc1'] = 'HLA-A*02:01'
data_formated['mhc2'] = 'B2M'

mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in snakemake.input.mhc_sequences])

mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Human', 'mhc_code'
].apply(tidytcells.mh.standardize, species='homosapiens')

mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
    mhc_sequences['species'] == 'Mouse', 'mhc_code'
].apply(tidytcells.mh.standardize, species='musmusculus')

data_formated['mhc_chain_1'] = data_formated['mhc1'].map(
    mhc_sequences.set_index('mhc_code')['canonical_sequence'].to_dict()
)

b2m_sequence = 'IQRTPKIQVYSRHPADIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM'
data_formated['mhc_chain_2'] = b2m_sequence

parser = PDBParser()
structure_1g4 = parser.get_structure('1g4', snakemake.input.tcr_1g4_structure)

alpha_chain_1g4 = ''.join(
    IUPACData.protein_letters_3to1[residue.get_resname().title()] for residue in structure_1g4[0]['D']
)
beta_chain_1g4 = ''.join(
    IUPACData.protein_letters_3to1[residue.get_resname().title()] for residue in structure_1g4[0]['E']
)

alpha_chain_c259 = alpha_chain_1g4[:93] + 'LY' + alpha_chain_1g4[95:]
beta_chain_c259 = beta_chain_1g4

cdr1_alpha, cdr2_alpha, cdr3_alpha = get_cdr_sequences(alpha_chain_c259)
cdr1_beta, cdr2_beta, cdr3_beta = get_cdr_sequences(beta_chain_c259)

data_formated['tcr_alpha_chain'] = alpha_chain_c259
data_formated['tcr_beta_chain'] = beta_chain_c259

data_formated['cdr1_alpha'] = cdr1_alpha
data_formated['cdr2_alpha'] = cdr2_alpha
data_formated['cdr3_alpha'] = cdr3_alpha

data_formated['cdr1_beta'] = cdr1_beta
data_formated['cdr2_beta'] = cdr2_beta
data_formated['cdr3_beta'] = cdr3_beta

data_formated.index = 'c259_' + data_formated.index

logger.info('Outputting data to %s', snakemake.output[0])
data_formated.to_csv(snakemake.output[0])
