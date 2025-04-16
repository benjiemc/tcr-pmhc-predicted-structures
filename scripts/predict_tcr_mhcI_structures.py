#!/usr/bin/env python
"""Run predictions of TCR:pMHC-I structures using TCRmodel2."""

import csv
import os
import subprocess

TCRMODEL2_ALPHAFOLD_DATA_DIR = os.environ['TCRMODEL2_ALPHAFOLD_DATA_DIR']

root_dir = os.getcwd()
os.chdir('third_party/tcrmodel2')

input_data = os.path.join(root_dir, 'data/interim/sampled_vdjdb_norm.csv')
output_dir = os.path.join(root_dir, 'data/interim/tcrmodel2_output/tcr-pmhc-Is')

with open(input_data, 'r') as fh:
    reader = csv.reader(fh)

    data = list(reader)

num_lines = len(data)
for i, (_, complex_id, alpha_chain_seq, beta_chain_seq, epitope_seq, mhc_seq) in enumerate(data, 1):
    print(f'STARTING COMPLEX: {complex_id} - {i} of {num_lines}', flush=True)

    if os.path.exists(os.path.join(output_dir, complex_id)):
        print(f'skipping {complex_id} as it has already been computed', flush=True)
        continue

    cmd = (
        'python '
        'run_tcrmodel2.py '
        f'--job_id={complex_id} '
        f'--output_dir={output_dir} '
        f'--tcra_seq={alpha_chain_seq} '
        f'--tcrb_seq={beta_chain_seq} '
        f'--pep_seq={epitope_seq} '
        f'--mhca_seq={mhc_seq} '
        '--relax_structures '
        f'--ori_db={TCRMODEL2_ALPHAFOLD_DATA_DIR}'
    )

    subprocess.run(cmd, shell=True)  # noqa: S602, PLW1510

    print('FINISHING COMPLEX: ', complex_id, flush=True)

print('WORKFLOW COMPLETED', flush=True)
