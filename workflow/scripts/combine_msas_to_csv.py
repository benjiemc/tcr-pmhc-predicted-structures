#!/usr/bin/env python
"""Combine a3m MSAs into one CSV formatted MSA."""

import re

import pandas as pd
from snakemake.script import snakemake

with open(snakemake.input.unpaired_msa, 'r') as fh:
    msa_unpaired = [line.strip() for line in fh if line]

msa_unpaired_grouped = [
    (msa_unpaired[i + 1], line.strip('>')) for i, line in enumerate(msa_unpaired) if line.startswith('>')
]
msa_unpaired_df = pd.DataFrame(msa_unpaired_grouped, columns=['sequence', 'comment'])
msa_unpaired_df['key'] = ''

with open(snakemake.input.paired_msa, 'r') as fh:
    msa_paired = [line.strip() for line in fh if line]

msa_paired_grouped = [(msa_paired[i + 1], line.strip('>')) for i, line in enumerate(msa_paired) if line.startswith('>')]
msa_paired_df = pd.DataFrame(msa_paired_grouped, columns=['sequence', 'comment'])

uniprot_entry_name_regex = re.compile(
    r'(?:tr|sp)\|(?:[A-Z0-9]{6,10})(?:_\d+)?\|(?:[A-Z0-9]{1,10}_)(?P<species>[A-Z0-9]{1,5})'
)
msa_paired_df['key'] = msa_paired_df['comment'].map(
    lambda comment: '' if (match := re.match(uniprot_entry_name_regex, comment)) is None else match.group('species')
)
msa_paired_df['key'] = msa_paired_df['key'].map(
    {species: str(i) if species else '' for i, species in enumerate(msa_paired_df['key'].unique())}
)

msa_df = pd.concat([msa_paired_df, msa_unpaired_df.iloc[1:]])
msa_df = msa_df[['sequence', 'key']]

msa_df.to_csv(snakemake.output[0], index=False)
