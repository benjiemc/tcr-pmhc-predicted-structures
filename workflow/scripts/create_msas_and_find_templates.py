#!/usr/bin/env python
"""Create MSAs and find templates for a sequence."""

import json
import os
import subprocess

from snakemake.script import snakemake

input_file_name = os.path.basename(snakemake.input[0])
input_dir_abs = os.path.join(os.getcwd(), os.path.dirname(snakemake.input[0]))
output_dir_abs = os.path.join(os.environ['TMPDIR'], snakemake.wildcards.name)

os.makedirs(output_dir_abs, exist_ok=True)

subprocess.run(  # noqa: S603
    [
        '/usr/bin/docker',
        'run',
        '-it',
        f'--volume={output_dir_abs}:/tmp/af_output',
        '--entrypoint',
        'bash',
        'alphafold3',
        '-c',
        f'chown -R {os.getuid()}:{os.getgid()} /tmp/af_output',
    ],
    check=True,
)

subprocess.run(  # noqa: S603
    [
        '/usr/bin/docker',
        'run',
        '-it',
        f'--user={os.getuid()}:{os.getgid()}',
        f'--volume={input_dir_abs}:/tmp/af_input',
        f'--volume={output_dir_abs}:/tmp/af_output',
        f'--volume={snakemake.params.alphafold_data_dir}:/public_databases',
        'alphafold3',
        'python',
        'run_alphafold.py',
        '--norun_inference',
        f'--json_path=/tmp/af_input/{input_file_name}',
        '--output_dir=/tmp/af_output',
    ],
    check=True,
)

os.makedirs(os.path.dirname(snakemake.output.unpaired_msa), exist_ok=True)
os.makedirs(os.path.dirname(snakemake.output.templates), exist_ok=True)

with open(os.path.join(output_dir_abs, snakemake.wildcards.name, f'{snakemake.wildcards.name}_data.json'), 'r') as fh:
    msa_data = json.load(fh)

    with open(snakemake.output.unpaired_msa, 'w') as out_fh:
        out_fh.write(msa_data['sequences'][0]['protein']['unpairedMsa'])

    with open(snakemake.output.paired_msa, 'w') as out_fh:
        out_fh.write(msa_data['sequences'][0]['protein']['pairedMsa'])

    os.mkdir(snakemake.output.templates)
    for i, template in enumerate(msa_data['sequences'][0]['protein']['templates']):
        with open(os.path.join(snakemake.output.templates, f'{i}.cif'), 'w') as out_fh:
            out_fh.write(template['mmcif'])

        with open(os.path.join(snakemake.output.templates, f'{i}_indices.json'), 'w') as out_fh:
            json.dump({key: template[key] for key in ('queryIndices', 'templateIndices')}, out_fh)
