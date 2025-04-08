import csv
import os
import subprocess

ALPHAFOLD_DATA_DIR = os.environ['ALPHAFOLD_DATA_DIR']

os.chdir("/home/imm/lina4096/code/tcrmodel2")

input_data = "/home/imm/lina4096/projects/tcr-pmhc-interface-analysis/data/sampled_vdjdb_norm.csv"
output_dir = "/home/imm/lina4096/projects/tcr-pmhc-interface-analysis/data/predicted-structures/tcr"

with open(input_data, 'r') as fh:
    reader = csv.reader(fh)
    header = next(reader) 
    data = [line for line in reader] 

num_lines = len(data)
for i, (_, complex_id, alpha_chain_seq, beta_chain_seq, _, _) in enumerate(data, 1):
    print(f"STARTING COMPLEX: {complex_id} - {i} of {num_lines}", flush=True)

    if os.path.exists(os.path.join(output_dir, complex_id)):
        print(f"skipping {complex_id} as it has already been computed", flush=True)
        continue

    cmd = ('python '
           'run_tcrmodel2_ub_tcr.py '
           f'--job_id={complex_id} '
           f'--output_dir={output_dir} '
           f'--tcra_seq={alpha_chain_seq} '
           f'--tcrb_seq={beta_chain_seq} '
           '--relax_structures '
           f'--ori_db={ALPHAFOLD_DATA_DIR}')

    subprocess.run(cmd, shell=True)
    
    print("FINISHING COMPLEX: ", complex_id, flush=True)

print("WORKFLOW COMPLETED", flush=True)
