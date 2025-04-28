"""Create AF3 input json from a CSV of TCR:pMHC sequence data."""

import argparse
import csv
import json
import logging
import os
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
parser.add_argument('--output', '-o', required=True, help='path to output json files')
parser.add_argument(
    '--mhc-sequences',
    nargs='+',
    required=True,
    help='paths to mhc sequence information to create pseudo sequences.',
)
parser.add_argument('--model-seeds', nargs='+', default=[1], help='')

add_logging_arguments(parser)

B2M_SEQUENCES = {
    'Human': 'MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPADIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM',
    'Mouse': (
        'MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEP'
        'KTVYWDRDM'
    ),
}
'''B2M sequences downloaded from UniProt.'''


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

    if not os.path.exists(args.output):
        logger.debug('Creating output dir: %s', args.output)
        os.mkdir(args.output)

    logger.info('Loading sequences from %s', args.input)
    with open(args.input, 'r') as fh:
        reader = csv.DictReader(fh)

        logger.info('Converting to JSON format and outputting to %s', args.output)
        for sequence in reader:
            logger.debug('Working on %s', sequence['id'])
            sequence['v_alpha'] = sequence['v_alpha'].split(';')[0]
            sequence['j_alpha'] = sequence['j_alpha'].split(';')[0]
            sequence['v_beta'] = sequence['v_beta'].split(';')[0]
            sequence['j_beta'] = sequence['j_beta'].split(';')[0]

            logger.debug('Stitching TCR alpha chain')
            tcr_alpha_seq = stitch_sequence(
                sequence['v_alpha'], sequence['j_alpha'], sequence['cdr3_alpha'], sequence['tcr_species']
            )

            if tcr_alpha_seq is None:
                logger.warning(
                    'Unable to stitch: %s-%s-%s, %s. Skipping...',
                    sequence['v_alpha'],
                    sequence['cdr3_alpha'],
                    sequence['j_alpha'],
                    sequence['tcr_species'],
                )
                continue

            logger.debug('Stitching TCR beta chain')
            tcr_beta_seq = stitch_sequence(
                sequence['v_beta'], sequence['j_beta'], sequence['cdr3_beta'], sequence['tcr_species']
            )

            if tcr_beta_seq is None:
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

            logger.debug('Looking up MHC chain 2 sequence')
            mhc_2_seq = (
                B2M_SEQUENCES[sequence['pmhc_species']]
                if sequence['mhc2'] == 'B2M'
                else mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc2']].iloc[0]['canonical_sequence']
            )

            output_name = os.path.join(args.output, f"{sequence['id']}.json")
            with open(output_name, 'w') as output_fh:
                logger.debug('Writing to %s', output_name)
                json.dump(
                    {
                        'name': sequence['id'],
                        'modelSeeds': args.model_seeds,
                        'sequences': [
                            {'protein': {'id': 'D', 'sequence': tcr_alpha_seq}},
                            {'protein': {'id': 'E', 'sequence': tcr_beta_seq}},
                            {'protein': {'id': 'C', 'sequence': sequence['peptide']}},
                            {'protein': {'id': 'A', 'sequence': mhc_1_seq}},
                            {'protein': {'id': 'B', 'sequence': mhc_2_seq}},
                        ],
                        'dialect': 'alphafold3',
                        'version': 1,
                    },
                    output_fh,
                )


if __name__ == '__main__':
    main()
