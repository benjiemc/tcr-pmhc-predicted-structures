#!/usr/bin/env python
"""Predict structure with Chai-1."""

import glob
import gzip
import hashlib
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path

from snakemake.script import snakemake

logging.basicConfig(level=logging.INFO)

with open(snakemake.input.sequences, 'r') as fh:
    lines = [line.strip() for line in fh if line]

sequences = {
    line.strip('>').split('|')[-1].replace('name=', ''): lines[i + 1]
    for i, line in enumerate(lines)
    if line.startswith('>')
}
sequence_ids = {
    name: {
        'id': (sequence_id := hashlib.sha256(sequence.encode()).hexdigest()),
        'id_short': sequence_id[:8],
        'sequence': sequence,
    }
    for name, sequence in sequences.items()
}

template_summaries = []
for template_summary_path in snakemake.input.templates_summaries:
    with open(template_summary_path, 'r') as fh:
        template_summaries.extend([line for line in fh if line])

new_template_summaries = []
for line in template_summaries:
    line_entries = line.split('\t')

    sequence_id = line_entries[0]
    chain_type = next(name for name, entry in sequence_ids.items() if entry['id'] == sequence_id)
    line_entries[0] = chain_type

    new_template_summaries.append('\t'.join(line_entries))

template_summaries = new_template_summaries

output_dir = os.path.dirname(snakemake.output.structures[0])
os.makedirs(output_dir, exist_ok=True)

with tempfile.TemporaryDirectory() as temp_dir:
    msa_dir = os.path.join(temp_dir, 'msas')
    os.mkdir(msa_dir)

    for msa in snakemake.input.msas:
        sequence_id_short = os.path.basename(msa).replace('.aligned.pqt', '')
        extended_name = next(entry['id'] for entry in sequence_ids.values() if entry['id_short'] == sequence_id_short)

        shutil.copy(msa, os.path.join(msa_dir, extended_name + '.aligned.pqt'))

    template_summary_path = os.path.join(temp_dir, 'templates_summary.m8')
    with open(template_summary_path, 'w') as fh:
        fh.write(''.join(template_summaries))

    template_dir = os.path.join(temp_dir, 'templates')
    os.mkdir(template_dir)
    os.environ['CHAI_TEMPLATE_CIF_FOLDER'] = template_dir

    for directory in snakemake.input.templates:
        name = os.path.basename(directory)

        for template_path in glob.glob(os.path.join(directory, '*.cif')):
            with open(template_path, 'rb') as fh_in:
                output_name = os.path.join(template_dir, name.upper() + '-' + os.path.basename(template_path) + '.gz')

                with gzip.open(output_name, 'wb') as fh_out:
                    shutil.copyfileobj(fh_in, fh_out)

    from chai_lab.chai1 import run_inference  # Import needs to happen after CHAI_TEMPLATE_CIF_FOLDER is set

    for seed in snakemake.config['seed']:
        temp_output_dir = os.path.join(temp_dir, f'output_{seed}')
        os.mkdir(temp_output_dir)

        candidates = run_inference(
            fasta_file=Path(snakemake.input.sequences),
            output_dir=Path(temp_output_dir),
            seed=int(seed),
            msa_directory=Path(msa_dir),
            template_hits_path=Path(template_summary_path),
        )

        for i, structure_path in enumerate(candidates.cif_paths):
            output_name = os.path.join(output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{i}_model.cif')
            shutil.copy(structure_path, output_name)

            summary_confidence_names = os.path.join(
                output_dir,
                f'{snakemake.wildcards.name}_seed-{seed}_sample-{i}_summary_confidences.json',
            )
            with open(summary_confidence_names, 'w') as fh:
                json.dump(
                    {
                        'aggregate_score': candidates.ranking_data[i].aggregate_score.cpu().numpy().tolist(),
                        'ptm': candidates.ranking_data[i].ptm_scores.complex_ptm.cpu().numpy().tolist(),
                        'iptm': candidates.ranking_data[i].ptm_scores.interface_ptm.cpu().numpy().tolist(),
                        'per_chain_ptm': candidates.ranking_data[i].ptm_scores.per_chain_ptm.cpu().numpy().tolist(),
                        'per_chain_pair_iptm': candidates.ranking_data[i]
                        .ptm_scores.per_chain_pair_iptm.cpu()
                        .numpy()
                        .tolist(),
                        'has_inter_chain_clashes': candidates.ranking_data[i]
                        .clash_scores.has_inter_chain_clashes.cpu()
                        .numpy()
                        .tolist(),
                        'chain_chain_clashes': candidates.ranking_data[i]
                        .clash_scores.chain_chain_clashes.cpu()
                        .numpy()
                        .tolist(),
                    },
                    fh,
                )

            confidence_name = os.path.join(
                output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{i}_confidences.json'
            )
            with open(confidence_name, 'w') as fh:
                json.dump(
                    {
                        'pae': candidates.pae[i].numpy().tolist(),
                        'pde': candidates.pde[i].numpy().tolist(),
                        'plddt': candidates.plddt[i].numpy().tolist(),
                    },
                    fh,
                )
