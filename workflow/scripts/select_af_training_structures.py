#!/usr/bin/env python
"""Select structures present in the AlphaFold training data and generate negatives."""

import logging
import os

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level='info')

ALPHAFOLD3_TRAINING_CUTOFF = '2021-09-30'


def get_sequences(path: str) -> dict[str, str]:
    """Get sequences for all chain in protein."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(
        os.path.basename(path).replace('.pdb', ''), os.path.join(snakemake.input.structures, path)
    )

    return {
        chain.id: ''.join(
            [IUPACData.protein_letters_3to1[res.get_resname().title()] for res in chain if res.id[0] == ' ']
        )
        for chain in structure[0]
    }


logger.info('Loading structure summary from %s', snakemake.input.structures_summary)
processed_structures = pd.read_csv(snakemake.input.structures_summary)

logger.info('Loading STCRDab summary from %s', snakemake.input.stcrdab)
stcrdab_summary = pd.read_csv(os.path.join(snakemake.input.stcrdab, 'db_summary.dat'), delimiter='\t')

# Known bug in STCRDab where entries that do get a date get assigned '01/08/09' (1909-01-08) as a date
stcrdab_summary = stcrdab_summary[stcrdab_summary['date'] != '01/08/09'].copy()
stcrdab_summary['date'] = pd.to_datetime(stcrdab_summary['date'], format=r'%m/%d/%y')

processed_structures = processed_structures.merge(
    stcrdab_summary[['pdb', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2', 'date']],
    how='left',
    on=['pdb', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2'],
)

logger.info('Selecting structures before %s', ALPHAFOLD3_TRAINING_CUTOFF)
in_distribution = processed_structures[processed_structures['date'] <= ALPHAFOLD3_TRAINING_CUTOFF]

logger.info('Selecting class I structures')
in_distribution = in_distribution.query("mhc_type == 'MH1'").copy()

logger.info('Extracting sequences from PDB files')
sequences = in_distribution['path'].apply(get_sequences)

filled_sequences = in_distribution.filter(regex=r'chain\d?(_x)?$').apply(
    lambda row: row.map(sequences[row.name].get), axis=1
)
filled_sequences = filled_sequences.dropna(axis='columns', how='all')
filled_sequences.columns = [col_name.replace('_x', '') + '_sequence' for col_name in filled_sequences.columns]

in_distribution = in_distribution.merge(filled_sequences, how='inner', left_index=True, right_index=True)

logger.info('Generating shuffled pairs')
rng = np.random.default_rng(snakemake.params.seed)

in_distribution['collated_cdrs'] = in_distribution[
    [
        'cdr1_alpha',
        'cdr2_alpha',
        'cdr3_alpha',
        'cdr1_beta',
        'cdr2_beta',
        'cdr3_beta',
    ]
].apply('-'.join, axis=1)
in_distribution['species'] = in_distribution['species'].fillna('Unknown')

peptides = in_distribution['peptide'].unique()

tcr_columns = [
    'pdb',
    'Achain',
    'Bchain',
    'cdr1_alpha',
    'cdr2_alpha',
    'cdr3_alpha',
    'cdr1_beta',
    'cdr2_beta',
    'cdr3_beta',
    'v_alpha',
    'v_beta',
    'species',
    'Achain_sequence',
    'Bchain_sequence',
]

pmhc_columns = [
    'pdb',
    'antigen_chain',
    'mhc_chain1',
    'mhc_chain2',
    'mhc_type',
    'peptide',
    'mhc1',
    'mhc2',
    'species',
    'antigen_chain_sequence',
    'mhc_chain1_sequence',
    'mhc_chain2_sequence',
]

negative_data = []

for peptide in peptides:
    in_group = in_distribution[in_distribution['peptide'] == peptide]
    out_group = in_distribution[
        (in_distribution['peptide'] != peptide) & (~in_distribution['collated_cdrs'].isin(in_group['collated_cdrs']))
    ]

    negatives = (out_group.sample(len(in_group), random_state=rng) if len(out_group) >= (len(in_group)) else out_group)[
        tcr_columns
    ]

    negatives = negatives.rename({'species': 'tcr_species', 'pdb': 'tcr_pdb'}, axis='columns')

    negatives = pd.concat(
        [
            negatives,
            pd.DataFrame(
                [in_group[pmhc_columns].iloc[0].to_numpy()] * len(negatives),
                columns=in_group[pmhc_columns].iloc[0].index,
                index=negatives.index,
            ),
        ],
        axis=1,
    )
    negatives['species'] = negatives[['tcr_species', 'species']].apply(
        lambda row: '/'.join({row.tcr_species, row.species}), axis=1
    )
    negatives = negatives.drop('tcr_species', axis='columns')

    negatives = negatives.rename({'pdb': 'pmhc_pdb'}, axis='columns')

    negative_data.append(negatives)

negative_data = pd.concat(negative_data).reset_index(drop=True)

negative_data['id'] = (
    negative_data['tcr_pdb']
    + '_'
    + negative_data['Achain']
    + negative_data['Bchain']
    + '_'
    + negative_data['pmhc_pdb']
    + '_'
    + negative_data['antigen_chain']
    + negative_data['mhc_chain1'].fillna('')
    + negative_data['mhc_chain2'].fillna('')
)
negative_data = negative_data.set_index('id')

negative_data['label'] = 0

logger.info('Merging selected structures and shuffled pairs')
in_distribution['id'] = (
    in_distribution['pdb']
    + '_'
    + in_distribution['Achain']
    + in_distribution['Bchain']
    + in_distribution['antigen_chain']
    + in_distribution['mhc_chain1'].fillna('')
    + in_distribution['mhc_chain2'].fillna('')
)
in_distribution = in_distribution.set_index('id')

in_distribution['label'] = 1

in_distribution = in_distribution.drop(['collated_cdrs', 'path', 'date'], axis='columns')

in_distribution['tcr_pdb'] = in_distribution['pdb']
in_distribution['pmhc_pdb'] = in_distribution['pdb']
in_distribution = in_distribution.drop('pdb', axis='columns')

output_data = pd.concat([in_distribution, negative_data])

logger.info('Writing output to %s', snakemake.output[0])
output_data.to_csv(snakemake.output[0])
