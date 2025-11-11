"""Extract TCR:pMHC relevant portion from the confidence metrics of a structure prediction."""

import argparse
import json
import logging
import sys

import numpy as np

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.confidence_scores import extract_confidence_score, get_sequence_information

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('confidence_path', help='paths to the confidence in numpy txt format')
parser.add_argument('sequence_path', help='path to the input sequence jsons')
parser.add_argument('--output', '-o', required=True, help='path to the output numpy array (txt format)')
parser.add_argument(
    '--score',
    choices=['atom_plddt', 'plddt', 'pae', 'pde'],
    required=True,
    help='score to extract from the input data',
)
parser.add_argument(
    '--binding-only',
    action='store_true',
    help=(
        'extract the confidence scores just isolated to the binding interface (CDRs, peptide, and MHC-TCR contact '
        'positions) rather than the full TCR and MHC chains'
    ),
)

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading sequences')
    with open(args.sequence_path, 'r') as fh:
        sequences = json.load(fh)

    sequences_flat, chain_ids, sequence_ids, insert_codes = get_sequence_information(sequences)

    logger.info('Loading confidence score')
    confidence = np.loadtxt(args.confidence_path)

    logger.info('Extracting TCR:pMHC confidence scores')
    selected_confidence = extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        confidence,
        score_type=args.score,
        mhc_type=sequences['mhc_type'],
    )

    logger.info('Outputting confidence to %s', args.output)
    np.savetxt(args.output, selected_confidence)


if __name__ == '__main__':
    main()
