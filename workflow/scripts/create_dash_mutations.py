#!/usr/bin/env python
"""Identify TCRs with motif from cluster and perform and alanine scan on the CDR3 regions."""

import logging

import anarci
import pandas as pd
import tidytcells
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_CDR3

logger = logging.getLogger()
setup_logger(logger, snakemake.config['log_level'])


def get_cdr3(sequence: str) -> str:
    """Get CDR3 from TCR chain sequence."""
    numbering, _ = anarci.number(sequence)
    cdr3 = ''.join(amino_acid for (seq_id, _), amino_acid in numbering if amino_acid != '-' and seq_id in IMGT_CDR3)
    return tidytcells.junction.standardise(cdr3)


def alanine_scan_cdr3(cdr3: str) -> list[str]:
    """Perform an alanine scane on all CDR3 residues."""
    new_cdr3s = []
    for i, residue in enumerate(cdr3[1:-1], 1):
        if residue == 'A':
            continue

        new_cdr3 = cdr3[:i] + 'A' + cdr3[i + 1 :]
        new_cdr3s.append(new_cdr3)

    return new_cdr3s


logger.info('Loading sequences from %s', snakemake.input.sequences)
dash_full_sequences = pd.read_csv(snakemake.input.sequences, index_col='id')

logger.info('Selecting sequences specific for %s', snakemake.params.peptide)
dash_full_sequences_epitope = dash_full_sequences[dash_full_sequences['peptide'] == snakemake.params.peptide].copy()

logger.info('Loading clusters from %s', snakemake.input.clusters)
with open(snakemake.input.clusters, 'r') as fh:
    clusters = dict(line.strip().split() for line in fh if line)

dash_full_sequences_epitope['cluster'] = dash_full_sequences_epitope.index.map(clusters)

logger.debug('Selecting only clustered data points')
dash_full_sequences_epitope_clustered = dash_full_sequences_epitope[
    dash_full_sequences_epitope['cluster'] != 'None'
].copy()

logger.info('Extract CDR3s from sequences')
dash_full_sequences_epitope_clustered['cdr3_alpha'] = dash_full_sequences_epitope_clustered['tcr_alpha_chain'].apply(
    get_cdr3
)
dash_full_sequences_epitope_clustered['cdr3_beta'] = dash_full_sequences_epitope_clustered['tcr_beta_chain'].apply(
    get_cdr3
)

logger.info('Creating mutations')
new_entries = []
for cluster in dash_full_sequences_epitope_clustered['cluster'].unique():
    logger.debug('Working on cluster %s', cluster)

    cdr3_alpha, cdr3_beta, _ = (
        dash_full_sequences_epitope_clustered[dash_full_sequences_epitope_clustered['cluster'] == cluster]
        .value_counts(['cdr3_alpha', 'cdr3_beta'])
        .reset_index()
        .iloc[0]
    )
    logger.debug('Top CDR3alpha: %s and CDR3beta: %s', cdr3_alpha, cdr3_beta)

    selected_tcr_pmhc = dash_full_sequences_epitope_clustered[
        (dash_full_sequences_epitope_clustered['cdr3_alpha'] == cdr3_alpha)
        & (dash_full_sequences_epitope_clustered['cdr3_beta'] == cdr3_beta)
    ].iloc[0]
    logger.debug('Selecting %s to mutate', selected_tcr_pmhc.name)

    new_chains = {}
    for chain_type in 'alpha', 'beta':
        logger.debug('Creating alanine substituitions for the %s-chain', chain_type)
        new_cdr3s = alanine_scan_cdr3(selected_tcr_pmhc[f'cdr3_{chain_type}'])

        logger.debug('Creating new full length sequnences')
        start_idx = selected_tcr_pmhc[f'tcr_{chain_type}_chain'].find(selected_tcr_pmhc[f'cdr3_{chain_type}'])
        end_idx = start_idx + len(selected_tcr_pmhc[f'cdr3_{chain_type}'])

        new_chains[chain_type] = [
            selected_tcr_pmhc[f'tcr_{chain_type}_chain'][:start_idx]
            + cdr3
            + selected_tcr_pmhc[f'tcr_{chain_type}_chain'][end_idx:]
            for cdr3 in new_cdr3s
        ]

    num = 1
    for chain_type, sequences in new_chains.items():
        for sequence in sequences:
            new_entry = {}
            for entry_name in (
                'tcr_alpha_chain',
                'tcr_beta_chain',
                'peptide',
                'mhc_type',
                'mhc_chain_1',
                'mhc_chain_2',
            ):
                new_entry[entry_name] = sequence if chain_type in entry_name else selected_tcr_pmhc[entry_name]

            new_entry['id'] = selected_tcr_pmhc.name + f'_mut_{num}'
            new_entries.append(new_entry)
            num += 1

new_entries = pd.DataFrame(new_entries)
new_entries = new_entries.set_index('id')

logger.info('Outputting mutated sequences to %s', snakemake.output[0])
new_entries.to_csv(snakemake.output[0])
