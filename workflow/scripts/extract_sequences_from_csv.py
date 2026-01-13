#!/usr/bin/env python
"""Script to extract sequences entries from CSVs to jsons."""

import json
import os

import pandas as pd
from snakemake.script import snakemake

df = pd.read_csv(snakemake.input[0])
output_df = df[
    [
        column_name
        for column_name in df.columns
        if column_name.endswith('_sequence') or column_name in ('id', 'mhc_type')
    ]
]

output_df.columns = [
    {
        'Achain_sequence': 'tcr_alpha_chain',
        'Bchain_sequence': 'tcr_beta_chain',
        'antigen_chain_sequence': 'peptide',
        'mhc_chain1_sequence': 'mhc_chain_1',
        'mhc_chain2_sequence': 'mhc_chain_2',
    }.get(column_name, column_name)
    for column_name in output_df.columns
]

output_dir = os.path.dirname(snakemake.output[0])
os.makedirs(output_dir, exist_ok=True)

for entry in [{k: v for k, v in row.items() if pd.notna(v)} for row in output_df.to_dict(orient='records')]:
    with open(os.path.join(output_dir, f"{entry['id']}.json"), 'w') as fh:
        json.dump(entry, fh)
