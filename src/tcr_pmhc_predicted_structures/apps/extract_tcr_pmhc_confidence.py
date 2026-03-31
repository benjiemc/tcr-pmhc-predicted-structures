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
    '--aggregate-atom-level-confidence',
    action='store_true',
    help='aggregate atom level confidence by averaging per residue',
)
parser.add_argument(
    '--portion',
    choices=['all', 'cropped', 'binding'],
    default='all',
    help=(
        "Extract the confidence scores for a portion of the complex. 'cropped' refers to the TCR-variable domain, "
        "peptide, and antigen binding domain. 'binding' is isolated to the binding interface (CDRs, peptide, and "
        "MHC-TCR contact positions) rather than the full TCR and MHC chains. (Default: 'all')"
    ),
)

padding_group = parser.add_argument_group('Padding')
padding_group.add_argument(
    '--pad',
    action='store_true',
    help='option to add padding to confidence scores being extract',
)
padding_group.add_argument(
    '--cdr1-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 1as to (Default: 8)',
)
padding_group.add_argument(
    '--cdr2-alpha-length',
    default=8,
    type=int,
    help='Length to pad CDR 2as to (Default: 8)',
)
padding_group.add_argument(
    '--cdr3-alpha-length',
    default=24,
    type=int,
    help='Length to pad CDR 3as to (Default: 24)',
)
padding_group.add_argument(
    '--cdr1-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 1bs to (Default: 8)',
)
padding_group.add_argument(
    '--cdr2-beta-length',
    default=8,
    type=int,
    help='Length to pad CDR 2bs to (Default: 8)',
)
padding_group.add_argument(
    '--cdr3-beta-length',
    default=24,
    type=int,
    help='Length to pad CDR 3bs to (Default: 24)',
)
padding_group.add_argument(
    '--peptide-length',
    default=12,
    type=int,
    help='Length to pad peptides to (Default: 12)',
)
padding_group.add_argument(
    '--pad-direction',
    choices=['centre', 'left', 'right'],
    default='centre',
    help="Direction to pad sequences (Default: 'centre')",
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

    if args.pad:
        logger.debug('Creating padding config')
        pad_config = {
            'direction': args.pad_direction,
            'cdr1_alpha_length': args.cdr1_alpha_length,
            'cdr2_alpha_length': args.cdr2_alpha_length,
            'cdr3_alpha_length': args.cdr3_alpha_length,
            'cdr1_beta_length': args.cdr1_beta_length,
            'cdr2_beta_length': args.cdr2_beta_length,
            'cdr3_beta_length': args.cdr3_beta_length,
            'peptide_length': args.peptide_length,
        }

    else:
        pad_config = None

    logger.info('Extracting TCR:pMHC confidence scores')
    selected_confidence = extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        confidence,
        portion=args.portion,
        score_type=args.score,
        mhc_type=sequences['mhc_type'],
        padding=pad_config,
        avg_atoms_for_residue=args.aggregate_atom_level_confidence,
    )

    logger.info('Outputting confidence to %s', args.output)
    np.savetxt(args.output, selected_confidence)


if __name__ == '__main__':
    main()
