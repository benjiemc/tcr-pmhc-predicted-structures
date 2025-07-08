"""Create full length sequences from TCR and MHC gene annotations."""

import argparse
import csv
import logging
import sys

import pandas as pd
import tidytcells

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.utils import (
    assign_species,
    mhc_slug_to_code,
    stitch_sequence,
)

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input', help='path to input csv')
parser.add_argument('--output', '-o', required=True, help='path to output csv')
parser.add_argument(
    '--mhc-sequences',
    nargs='+',
    required=True,
    help='paths to mhc sequence information to create pseudo sequences.',
)

add_logging_arguments(parser)

B2M_SEQUENCES = {
    'Human': 'MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPADIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM',
    'Mouse': (
        'MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEP'
        'KTVYWDRDM'
    ),
}
'''B2M sequences downloaded from UniProt.'''


def create_full_length_sequences(sequences: list[dict], mhc_sequences: pd.DataFrame) -> list[dict]:
    """Create full length TCR:pMHC sequences from V and J genes, CDR3s, MHC allele, and peptide sequence."""
    sequences_full = []

    for sequence in sequences:
        logger.debug('Working on %s', sequence['index'])
        sequence['v_alpha'] = sequence['v_alpha'].split(';')[0]
        sequence['j_alpha'] = sequence['j_alpha'].split(';')[0]
        sequence['v_beta'] = sequence['v_beta'].split(';')[0]
        sequence['j_beta'] = sequence['j_beta'].split(';')[0]

        logger.debug('Stitching TCR alpha chain')
        for species in sequence['species'].split('/'):
            tcr_alpha_seq = stitch_sequence(sequence['v_alpha'], sequence['j_alpha'], sequence['cdr3_alpha'], species)

            if tcr_alpha_seq is not None:
                break

        else:
            logger.warning(
                'Unable to stitch: %s-%s-%s, %s. Skipping...',
                sequence['v_alpha'],
                sequence['cdr3_alpha'],
                sequence['j_alpha'],
                sequence['tcr_species'],
            )
            continue

        logger.debug('Stitching TCR beta chain')
        for species in sequence['species'].split('/'):
            tcr_beta_seq = stitch_sequence(sequence['v_beta'], sequence['j_beta'], sequence['cdr3_beta'], species)

            if tcr_beta_seq is not None:
                break

        else:
            logger.warning(
                'Unable to stitch: %s-%s-%s, %s. Skipping...',
                sequence['v_beta'],
                sequence['cdr3_beta'],
                sequence['j_beta'],
                sequence['tcr_species'],
            )
            continue

        logger.debug('Looking up MHC chain 1 sequence')
        mhc_1_seq = mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc1']].iloc[0]['canonical_sequence']
        mhc_species = mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc1']].iloc[0]['species']

        logger.debug('Looking up MHC chain 2 sequence')
        mhc_2_seq = (
            B2M_SEQUENCES[mhc_species]
            if sequence['mhc2'] == 'B2M'
            else mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc2']].iloc[0]['canonical_sequence']
        )

        sequences_full.append(
            {
                'id': sequence['index'],
                'tcr_alpha_chain': tcr_alpha_seq,
                'tcr_beta_chain': tcr_beta_seq,
                'peptide': sequence['peptide'],
                'mhc_chain_1': mhc_1_seq,
                'mhc_chain_2': mhc_2_seq,
                'mhc_type': sequence['mhc_type'],
            }
        )

    return sequences_full


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading MHC canonical sequences')
    mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in args.mhc_sequences])

    mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
    mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

    mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Human', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='homosapiens')

    mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Mouse', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='musmusculus')

    logger.info('Loading sequences from %s', args.input)
    with open(args.input, 'r') as fh:
        sequences = list(csv.DictReader(fh))

    logger.info('Extending sequences')
    full_sequences = create_full_length_sequences(sequences, mhc_sequences)

    logger.info('Writing sequences to %s', args.output)
    with open(args.output, 'w') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                'id',
                'tcr_alpha_chain',
                'tcr_beta_chain',
                'peptide',
                'mhc_chain_1',
                'mhc_chain_2',
                'mhc_type',
            ],
        )
        writer.writeheader()
        writer.writerows(full_sequences)


if __name__ == '__main__':
    main()
