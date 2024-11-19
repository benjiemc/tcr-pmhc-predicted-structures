import argparse

import anarci
import pandas as pd

IMGT_CDR1 = set(range(27, 38 + 1))
IMGT_CDR2 = set(range(56, 65 + 1))
IMGT_CDR3 = set(range(105, 117 + 1))

parser = argparse.ArgumentParser()

parser.add_argument('input_file')
parser.add_argument('--output', '-o', help='path to output csv')


def main():
    args = parser.parse_args()

    vdjdb = pd.read_csv(args.input_file)
    
    pseudo_seqs = []
    for i, row in vdjdb.iterrows():
        print('Working on:', row.complex_id, flush=True)
        pseudo_seq = {}

        for chain in 'alpha_chain_sequence', 'beta_chain_sequence':
            sequence = row[chain]
            numbered_sequence, _ = anarci.number(sequence)

            cdr1 = ''
            cdr2 = ''
            cdr3 = ''

            if not numbered_sequence:
                print('Anarci numbering failed...', flush=True)

            else:
                for (seq_id, _), res_name in numbered_sequence:
                    if res_name == '-':
                        continue

                    elif seq_id in IMGT_CDR1:
                        cdr1 += res_name

                    elif seq_id in IMGT_CDR2:
                        cdr2 += res_name

                    elif seq_id in IMGT_CDR3:
                        cdr3 += res_name

            pseudo_seq[chain] = '-'.join([cdr1, cdr2, cdr3])
        
        pseudo_seqs.append('-'.join([pseudo_seq['alpha_chain_sequence'], pseudo_seq['beta_chain_sequence']]))

    vdjdb['tcr_pseudo_sequence'] = pseudo_seqs
    vdjdb.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
