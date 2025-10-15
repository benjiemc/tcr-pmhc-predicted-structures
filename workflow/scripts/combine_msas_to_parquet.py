#!/usr/bin/env python
"""Combine a3m MSAs into one parquet formatted MSA."""

import re

import pandas as pd
from snakemake.script import snakemake

CHAI_ALLOWED_DATABASES = {'uniref90', 'mgnify', 'bfd_uniclust', 'uniprot', 'query'}

with open(snakemake.input.unpaired_msa, 'r') as fh:
    msa_unpaired = [line.strip() for line in fh if line]

msa_unpaired_grouped = [
    (msa_unpaired[i + 1], line.strip('>')) for i, line in enumerate(msa_unpaired) if line.startswith('>')
]
msa_unpaired_df = pd.DataFrame(msa_unpaired_grouped, columns=['sequence', 'comment'])


msa_unpaired_df['source_database'] = msa_unpaired_df['comment'].map(lambda comment: comment.split()[0].split('_')[0])
msa_unpaired_df['source_database'] = msa_unpaired_df['source_database'].map(
    lambda database_name: 'uniprot' if re.fullmatch(r'[A-Z0-9]+', database_name) else database_name
)
msa_unpaired_df['source_database'] = msa_unpaired_df['source_database'].str.lower()
msa_unpaired_df.loc[~msa_unpaired_df['source_database'].isin(CHAI_ALLOWED_DATABASES), 'source_database'] = (
    'uniprot'  # Hack to avoid errors, might not be the correct DB
)

msa_unpaired_df['pairing_key'] = ''

with open(snakemake.input.paired_msa, 'r') as fh:
    msa_paired = [line.strip() for i, line in enumerate(fh) if line]

msa_paired_grouped = [(msa_paired[i + 1], line.strip('>')) for i, line in enumerate(msa_paired) if line.startswith('>')]
msa_paired_df = pd.DataFrame(msa_paired_grouped, columns=['sequence', 'comment'])

msa_paired_df['source_database'] = 'uniprot'
msa_paired_df.loc[msa_paired_df['comment'] == 'query', 'source_database'] = 'query'
msa_paired_df['source_database'] = msa_paired_df['source_database'].str.lower()

uniprot_entry_name_regex = re.compile(
    r'(?:tr|sp)\|(?:[A-Z0-9]{6,10})(?:_\d+)?\|(?:[A-Z0-9]{1,10}_)(?P<species>[A-Z0-9]{1,5})'
)
msa_paired_df['pairing_key'] = msa_paired_df['comment'].map(
    lambda comment: ('' if (match := re.match(uniprot_entry_name_regex, comment)) is None else match.group('species'))
)

msa_df = pd.concat([msa_paired_df, msa_unpaired_df[msa_unpaired_df['source_database'] != 'query']])
msa_df = msa_df[['sequence', 'source_database', 'pairing_key', 'comment']]

msa_df.to_parquet(snakemake.output[0], index=False)
