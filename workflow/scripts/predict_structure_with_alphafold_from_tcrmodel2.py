#!/usr/bin/env python
"""Predict structure with Alphafold from TCRmodel2."""

import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBIO, MMCIFParser
from snakemake.script import snakemake

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

with open(snakemake.input.sequences, 'r') as fh:
    lines = [line.strip() for line in fh if line]

sequences = [
    (line.strip('>').split('|')[-1].replace('name=', ''), lines[i + 1])
    for i, line in enumerate(lines)
    if line.startswith('>')
]

sequence_order = [sequence_type for sequence_type, _ in sequences]
sequence_ids = {
    sequence_type: {'id': hashlib.sha256(sequence.encode()).hexdigest()[:8], 'sequence': sequence}
    for sequence_type, sequence in sequences
}

output_dir = os.path.dirname(snakemake.output.structures_unrelaxed[0])

with tempfile.TemporaryDirectory() as tmp_dir:
    # Fasta sequence
    shutil.copy(snakemake.input.sequences, os.path.join(tmp_dir, 'sequences.fasta'))

    # Templates
    temp_template_dir = os.path.join(tmp_dir, 'templates')
    os.mkdir(temp_template_dir)

    alignmnet_files = []
    for sequence_type in sequence_order:
        if sequence_type == 'peptide':
            alignmnet_files.append('')
            continue

        template_dir = next(
            path for path in snakemake.input.templates if sequence_ids[sequence_type]['id'] in os.path.basename(path)
        )
        name = os.path.basename(template_dir)

        template_entries = defaultdict(dict)
        indices_paths = glob.glob(os.path.join(template_dir, '*_indices.json'))
        structure_paths = glob.glob(os.path.join(template_dir, '*.cif'))

        for indice_path in indices_paths:
            entry_name = re.search(r'(\d+)_indices.json', os.path.basename(indice_path)).group(1)
            template_entries[entry_name]['indices_path'] = indice_path

        for structure_path in structure_paths:
            entry_name = re.search(r'(\d+).cif', os.path.basename(structure_path)).group(1)
            template_entries[entry_name]['structure_path'] = structure_path

        template_entries = sorted(template_entries.items())

        with open(os.path.join(temp_template_dir, name + '_alignment.tsv'), 'w') as fh_out:
            fh_out.write('template_pdbfile\ttarget_to_template_alignstring\ttarget_len\ttemplate_len\tidentities\n')

            for template_name, template_info in template_entries:
                output_name = f'{name}_{template_name}.pdb'
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure('', template_info['structure_path'])

                chain = next(structure[0].get_chains())
                template_length = len(chain)
                valid_indices = {idx for idx, residue in enumerate(chain) if residue.id[0] == ' '}

                index_mapping = {}
                valid_indices = set()
                new_idx = 0
                for idx, residue in enumerate(chain):
                    if residue.id[0] != ' ':
                        continue

                    index_mapping[idx] = new_idx
                    valid_indices.add(idx)
                    new_idx += 1

                io = PDBIO()
                io.set_structure(structure)
                io.save(os.path.join(temp_template_dir, output_name))

                with open(template_info['indices_path'], 'r') as fh:
                    indices = json.load(fh)

                target_to_template_align_string = ';'.join(
                    [
                        f'{query_idx}:{index_mapping[template_idx]}'
                        for query_idx, template_idx in zip(
                            indices['queryIndices'], indices['templateIndices'], strict=True
                        )
                        if template_idx in valid_indices
                    ]
                )

                target_length = len(sequence_ids[sequence_type]['sequence'])

                fh_out.write(
                    f'/root/input/templates/{output_name}\t'
                    f'{target_to_template_align_string}\t'
                    f'{target_length}\t'
                    f'{template_length}\t'
                    '0\n'
                )

        alignmnet_files.append(f'/root/input/templates/{name}_alignment.tsv')

    alignmnet_files = ','.join(alignmnet_files)

    for seed in snakemake.config['seed']:
        temp_output_dir = os.path.join(tmp_dir, f'output_{seed}')
        os.mkdir(temp_output_dir)

        # MSAs
        msa_dir = os.path.join(temp_output_dir, 'sequences', 'msas')
        os.makedirs(msa_dir)

        for chain_id, sequence_type in zip(PDB_CHAIN_IDS, sequence_order, strict=False):
            msa_chain_dir = os.path.join(msa_dir, chain_id)
            os.mkdir(msa_chain_dir)

            if sequence_type == 'peptide':
                with open(os.path.join(msa_chain_dir, 'custom_MSA.sto'), 'w') as fh:
                    fh.write(f"# STOCKHOLM 1.0\nquery {sequence_ids['peptide']['sequence']}\n//\n")

            else:
                msa_path = next(
                    path for path in snakemake.input.msas if sequence_ids[sequence_type]['id'] in os.path.basename(path)
                )
                with open(msa_path, 'r') as fh_in, open(os.path.join(msa_chain_dir, 'custom_MSA.sto'), 'w') as fh_out:
                    fh_out.write('# STOCKHOLM 1.0\n')

                    while True:
                        line = fh_in.readline()
                        if not line:
                            break

                        if not line.startswith('>'):
                            continue

                        description = '_'.join(line.strip().strip('>').split())
                        content = fh_in.readline().strip().upper()

                        fh_out.write(f'{description} {content}\n')

                    fh_out.write('//\n')

        subprocess.run(  # noqa: S603
            [
                '/usr/bin/apptainer',
                'exec',
                '--nv',
                f'--bind={snakemake.params.tcrmodel2_alphafold_data_dir}:/root/data_dir',
                f'--bind={temp_output_dir}:/root/output',
                f'--bind={tmp_dir}:/root/input',
                'containers/tcrmodel2.sif',
                'run_alphafold_tcrmodel2.3',
                '--model_preset=multimer',
                f'--random_seed={seed}',
                '--num_multimer_predictions_per_model=1',
                '--models_to_relax=all',
                '--use_gpu_relax',
                '--use_custom_MSA_database=only',
                '--MSA_database='  # Needed not to throw an error
                '/data/databases/uniref90.tcrmhc.fasta,'
                '/data/databases/uniref90.tcrmhc.fasta,'
                '/data/databases/uniref90.tcrmhc.fasta,'
                '/data/databases/uniref90.tcrmhc.fasta,'
                '/data/databases/uniref90.tcrmhc.fasta',
                '--db_preset=reduced_dbs',  # Needed not to throw an error
                '--uniref90_database_path=/data/databases/uniref90.tcrmhc.fasta',  # Needed not to throw an error
                '--uniprot_database_path=/data/databases/uniprot.tcrmhc.fasta',  # Needed not to throw an error
                '--mgnify_database_path=/data/databases/mgnify.fasta',  # Needed not to throw an error
                '--small_bfd_database_path=/data/databases/small_bfd.tcrmhc.fasta',  # Needed not to throw an error
                '--pdb_seqres_database_path=/data/databases/pdb_seqres.txt',
                '--template_mmcif_dir=/root/data_dir/pdb_mmcif/mmcif_files/',
                '--obsolete_pdbs_path=/root/data_dir/pdb_mmcif/obsolete.dat',
                '--max_template_date=2100-01-01',
                '--use_custom_templates',
                f'--template_alignfile={alignmnet_files}',
                '--data_dir=/root/data_dir',
                '--fasta_paths=/root/input/sequences.fasta',
                '--output_dir=/root/output',
            ],
            check=True,
        )

        for sample in range(1, 5 + 1):
            shutil.copy(
                os.path.join(temp_output_dir, 'sequences', f'unrelaxed_model_{sample}_multimer_v3_pred_0.pdb'),
                os.path.join(output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample - 1}_unrelaxed.pdb'),
            )
            shutil.copy(
                os.path.join(temp_output_dir, 'sequences', f'relaxed_model_{sample}_multimer_v3_pred_0.pdb'),
                os.path.join(output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample - 1}_relaxed.pdb'),
            )

            confidence_scores = np.load(
                os.path.join(temp_output_dir, 'sequences', f'result_model_{sample}_multimer_v3_pred_0.pkl'),
                allow_pickle=True,
            )
            with open(
                os.path.join(
                    output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample - 1}_summary_confidences.json'
                ),
                'w',
            ) as fh:
                json.dump(
                    {
                        'ptm': float(confidence_scores['ptm']),
                        'iptm': float(confidence_scores['iptm']),
                        'ranking_confidence': float(confidence_scores['ranking_confidence']),
                        'max_pae': float(confidence_scores['max_predicted_aligned_error']),
                    },
                    fh,
                )

            with open(
                os.path.join(
                    output_dir, f'{snakemake.wildcards.name}_seed-{seed}_sample-{sample - 1}_confidences.json'
                ),
                'w',
            ) as fh:
                json.dump(
                    {
                        'pae': confidence_scores['predicted_aligned_error'].tolist(),
                        'plddt': confidence_scores['plddt'].tolist(),
                    },
                    fh,
                )
