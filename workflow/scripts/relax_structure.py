#!/usr/bin/env python
"""Relax strucuture using Amber provided by AlphaFold2/TCRmodel2."""

import os
import shutil
import subprocess
import tempfile

from snakemake.script import snakemake

with tempfile.TemporaryDirectory() as tmp_dir:
    input_dir = os.path.join(tmp_dir, 'input')
    os.mkdir(input_dir)

    output_dir = os.path.join(tmp_dir, 'output')
    os.mkdir(output_dir)

    input_name = os.path.basename(snakemake.input[0])
    shutil.copy(snakemake.input[0], os.path.join(input_dir, input_name))

    subprocess.run(  # noqa: S603
        [
            '/usr/bin/apptainer',
            'exec',
            '--nv',
            f'--bind={input_dir}:/root/input',
            f'--bind={output_dir}:/root/output',
            'containers/tcrmodel2.sif',
            'relax_structure',
            '--use_gpu_relax',
            f'--input_pdb=/root/input/{input_name}',
            '--output_dir=/root/output',
        ],
        check=True,
    )

    shutil.copy(os.path.join(output_dir, f'relaxed_{input_name}'), snakemake.output[0])
