#!/usr/bin/env python
"""Predict structures with Boltz."""

import glob
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np
from snakemake.script import snakemake

with tempfile.TemporaryDirectory() as tmp_dir:
    input_dir_abs = os.path.join(tmp_dir, snakemake.wildcards.name)
    os.mkdir(input_dir_abs)

    output_dir = os.path.join(os.getcwd(), os.path.dirname(snakemake.output.structures[0]))
    os.makedirs(output_dir, exist_ok=True)

    shutil.copy(snakemake.input.sequence, input_dir_abs)

    for msa_path in snakemake.input.msas:
        shutil.copy(msa_path, input_dir_abs)

    for template_dir in snakemake.input.templates:
        directory_name = os.path.basename(template_dir)
        new_directory_path = os.path.join(input_dir_abs, directory_name)
        os.mkdir(new_directory_path)

        for template in glob.glob(os.path.join(template_dir, '*.cif')):
            template_name = os.path.basename(template)
            shutil.copy(template, os.path.join(new_directory_path, f'{directory_name}_{template_name}'))

    subprocess.Popen(['/usr/bin/chmod', '+x', '-R', input_dir_abs])  # noqa: S603

    os.chdir(input_dir_abs)

    input_file_name = os.path.basename(snakemake.input.sequence)

    for seed in snakemake.config['seed']:
        temp_output_dir_abs = os.path.join(tmp_dir, snakemake.wildcards.name, f'output_{seed}')
        os.mkdir(temp_output_dir_abs)

        subprocess.run(  # noqa: S603
            [  # noqa: S607 allow boltz to be called without full path
                'boltz',
                'predict',
                f"--cache={os.path.join(os.environ['TMPDIR'], 'boltz_cache')}",
                f'--seed={seed}',
                '--diffusion_samples=5',
                '--max_parallel_samples=1',
                '--write_full_pae',
                '--write_full_pde',
                f'--out_dir={temp_output_dir_abs}',
                os.path.join(input_dir_abs, input_file_name),
            ],
            check=True,
        )

        for sample_num in range(5):
            prediction_dir = os.path.join(
                temp_output_dir_abs,
                f'boltz_results_{snakemake.wildcards.name}',
                'predictions',
                snakemake.wildcards.name,
            )

            shutil.copy(
                os.path.join(prediction_dir, f'{snakemake.wildcards.name}_model_{sample_num}.cif'),
                os.path.join(output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample_num}_model.cif'),
            )
            shutil.copy(
                os.path.join(prediction_dir, f'confidence_{snakemake.wildcards.name}_model_{sample_num}.json'),
                os.path.join(
                    output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample_num}_summary_confidences.json'
                ),
            )

            with open(
                os.path.join(
                    output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample_num}_confidences.json'
                ),
                'w',
            ) as fh:
                json.dump(
                    {
                        'pae': np.load(
                            os.path.join(prediction_dir, f'pae_{snakemake.wildcards.name}_model_{sample_num}.npz')
                        )['pae'].tolist(),
                        'pde': np.load(
                            os.path.join(prediction_dir, f'pde_{snakemake.wildcards.name}_model_{sample_num}.npz')
                        )['pde'].tolist(),
                        'plddt': np.load(
                            os.path.join(prediction_dir, f'plddt_{snakemake.wildcards.name}_model_{sample_num}.npz')
                        )['plddt'].tolist(),
                    },
                    fh,
                )
