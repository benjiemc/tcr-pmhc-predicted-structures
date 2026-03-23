#!/usr/bin/env python
"""Create MSA inputs from secondary datasets (ignoring entries that are already present in the primary ones)."""

import csv
import hashlib
import json
import os

from snakemake.script import snakemake

output_dir = os.path.dirname(snakemake.output[0])
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

past_ids = set()
with open(snakemake.input.sequence_key, 'r') as fh:
    past_ids = past_ids | {line.strip().split()[0] for line in fh if line}

if hasattr(snakemake.input, 'higher_priority_datasets'):
    for dataset_path in snakemake.input.higher_priority_datasets:
        with open(dataset_path, 'r') as fh:
            for entry in csv.DictReader(fh):
                for chain_name in 'tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2':
                    if chain_name in entry:
                        past_ids.add(hashlib.sha256(entry[chain_name].encode()).hexdigest()[:8])

with open(snakemake.input.sequences, 'r') as fh:
    reader = csv.DictReader(fh)
    for entry in reader:
        for sequence_name in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2'):
            if sequence_name not in entry:
                continue

            sequence = entry[sequence_name]
            if not sequence:
                continue

            sequence_id = hashlib.sha256(sequence.encode()).hexdigest()[:8]
            if sequence_id in past_ids:
                continue

            with open(os.path.join(output_dir, sequence_id + '.json'), 'w') as output_fh:
                json.dump(
                    {
                        'name': sequence_id,
                        'modelSeeds': [1],
                        'sequences': [
                            {'protein': {'id': 'A', 'sequence': sequence}},
                        ],
                        'dialect': 'alphafold3',
                        'version': 1,
                    },
                    output_fh,
                )

            past_ids.add(sequence_id)
