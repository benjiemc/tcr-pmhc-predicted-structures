#!/usr/bin/env python
"""Create MSAs and find templates for a sequence."""

import json
import os
import subprocess
import tempfile

from snakemake.script import snakemake

NUM_PARALLEL_AF_PROCESSES = 4
'''AF3 runs 4 processes at once for building the MSAs.'''

with tempfile.TemporaryDirectory() as temp_dir:
    subprocess.run(  # noqa: S603
        [
            '/usr/bin/apptainer',
            'exec',
            f'--bind={snakemake.params.alphafold_data_dir}:/root/public_databases',
            (
                f"--bind={snakemake.params.pdb_database_path.replace('${DB_DIR}', snakemake.params.alphafold_data_dir)}"
                ":/root/pdb"
            ),
            'containers/alphafold3.sif',
            'python',
            'third_party/alphafold3/run_alphafold.py',
            '--norun_inference',
            f'--jackhmmer_n_cpu={max(snakemake.threads // NUM_PARALLEL_AF_PROCESSES, 1)}',
            f'--nhmmer_n_cpu={max(snakemake.threads // NUM_PARALLEL_AF_PROCESSES, 1)}',
            '--db_dir=/root/public_databases',
            '--pdb_database_path=/root/pdb',
            f'--uniref90_database_path={snakemake.params.uniref90_database_path}',
            f'--mgnify_database_path={snakemake.params.mgnify_database_path}',
            f'--small_bfd_database_path={snakemake.params.small_bfd_database_path}',
            f'--uniprot_cluster_annot_database_path={snakemake.params.uniprot_cluster_annot_database_path}',
            f'--seqres_database_path={snakemake.params.seqres_database_path}',
            f'--json_path={snakemake.input[0]}',
            f'--output_dir={temp_dir}',
        ],
        check=True,
    )

    os.makedirs(os.path.dirname(snakemake.output.unpaired_msa), exist_ok=True)
    os.makedirs(os.path.dirname(snakemake.output.templates), exist_ok=True)

    with open(os.path.join(temp_dir, snakemake.wildcards.name, f'{snakemake.wildcards.name}_data.json'), 'r') as fh:
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
