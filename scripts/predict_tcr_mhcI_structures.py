import csv
import os
import subprocess

ALPHAFOLD_DATA_DIR = os.environ['ALPHAFOLD_DATA_DIR']

os.chdir("/home/benjiemc/code/tcrmodel2")

input_data = "/home/benjiemc/projects/tcr-pmhc-predicted-structures/scripts/input.txt"
output_dir = "/home/benjiemc/projects/tcr-pmhc-predicted-structures/data/predicted-structures"

with open(input_data, 'r') as fh:
    reader = csv.reader(fh)
    
    data = [line for line in reader] 

num_lines = len(data)
for i, (_, complex_id, alpha_chain_seq, beta_chain_seq, epitope_seq, mhc_seq) in enumerate(data, 1):
    print(f"STARTING COMPLEX: {complex_id} - {i} of {num_lines}", flush=True)

    if os.path.exists(os.path.join(output_dir, complex_id)):
        print(f"skipping {complex_id} as it has already been computed", flush=True)
        continue

    cmd = ('python '
           'run_tcrmodel2.py '
           f'--job_id={complex_id} '
           f'--output_dir={output_dir} '
           f'--tcra_seq={alpha_chain_seq} '
           f'--tcrb_seq={beta_chain_seq} '
           f'--pep_seq={epitope_seq} '
           f'--mhca_seq={mhc_seq} '
           '--relax_structures '
           f'--ori_db={ALPHAFOLD_DATA_DIR}')

    subprocess.run(cmd, shell=True)
    
    print("FINISHING COMPLEX: ", complex_id, flush=True)

print("WORKFLOW COMPLETED", flush=True)
