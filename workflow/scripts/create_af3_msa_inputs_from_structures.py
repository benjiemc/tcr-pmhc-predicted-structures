#!/usr/bin/env python
"""Create AF3 MSA inputs from CSVs of sequences taken from structure datasets."""

import csv
import hashlib
import json
import os

from snakemake.script import snakemake

output_dir = os.path.dirname(snakemake.output[0])
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

past_ids = set()

if hasattr(snakemake.input, 'sequence_key'):
    with open(snakemake.input.sequence_key, 'r') as fh:
        past_ids = past_ids | {line.strip().split()[0] for line in fh if line}

if hasattr(snakemake.input, 'additional_structures'):
    with open(snakemake.input.additional_structures, 'r') as fh:
        for entry in csv.DictReader(fh):
            for sequence_name in ('Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence'):
                if sequence_name not in entry:
                    continue

                sequence = entry[sequence_name]
                if not sequence:
                    continue

                past_ids.add(hashlib.sha256(sequence.encode()).hexdigest()[:8])

with open(snakemake.input.structures, 'r') as fh:
    for entry in csv.DictReader(fh):
        for sequence_name in ('Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence'):
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
