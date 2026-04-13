#! /usr/bin/env python
"""Aggregate and normalise contact maps from PDB structures."""

import logging
import os

import pandas as pd
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Aggregating contacts')
contacts = []
for path in snakemake.input:
    logger.debug('Loading contacts from %s', path)
    contact = pd.read_csv(path, dtype={'resi_pmhc': str})

    name = os.path.basename(path).replace('_contacts.csv', '')
    contact['name'] = name

    contacts.append(contact)

contacts = pd.concat(contacts)

logger.info('Normalising contacts')
contacts = contacts.drop_duplicates(
    [
        'name',
        'chain_type_tcr',
        'residue_seq_id_tcr',
        'residue_insert_code_tcr',
        'chain_type_pmhc',
        'residue_seq_id_pmhc',
        'residue_insert_code_pmhc',
    ]
)
contacts['cdr_type'] = (
    'CDR'
    + contacts['cdr_tcr'].apply(int).apply(str)
    + contacts['chain_type_tcr'].map({'Achain': 'alpha', 'Bchain': 'beta'})
)

norm_columns = [
    'cdr_type',
    'relative_pos_centre_tcr',
    'mhc_type',
    'chain_type_pmhc',
    'resi_pmhc',
    'relative_pos_centre_pmhc',
]
contacts = contacts[contacts['chain_type_pmhc'].isin(['mhc_chain1', 'mhc_chain2', 'antigen_chain'])]

contact_counts = contacts.value_counts(norm_columns, dropna=False, normalize=True)

logger.info('Outputting results to %s', snakemake.output[0])
contact_counts = contact_counts.to_frame().reset_index()
contact_counts = contact_counts.sort_values(norm_columns).reset_index(drop=True)

contact_counts.to_csv(snakemake.output[0], index=False)
