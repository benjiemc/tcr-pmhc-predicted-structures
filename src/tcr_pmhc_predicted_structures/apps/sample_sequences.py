"""Sample sequence data for structure predictions."""

import argparse
import logging
import sys

import numpy as np
import pandas as pd

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input', help='path to input csv')
parser.add_argument('--output', '-o', required=True, help='path to output csv')
parser.add_argument('--sample-size', type=int, default=1000, help='size of the sample to take (Default: 1000)')
parser.add_argument('--seed', type=int, default=None, help='seed for random sampling')

background_group = parser.add_argument_group('Background')
background_group.add_argument(
    '--generate-background',
    action='store_true',
    help='option to create a background by randomly mixing TCRs and pMHCs',
)
background_group.add_argument(
    '--background-proportion',
    default=1,
    type=int,
    help=(
        'ratio of true pairings to random pairings for the background (Default: 1, meaning equal true pairings to '
        'background)'
    ),
)

add_logging_arguments(parser)


def generate_background(
    sequences: pd.DataFrame, background_ratio: float = 1.0, random_state: np.random.BitGenerator | None = None
) -> pd.DataFrame:
    sequences = sequences.copy()
    sequences['collated_cdrs'] = sequences.filter(regex='^cdr[1-3]_(alpha|beta)').apply('-'.join, axis=1)

    tcr_columns = [
        'id',
        'cdr1_alpha',
        'cdr2_alpha',
        'cdr3_alpha',
        'cdr1_beta',
        'cdr2_beta',
        'cdr3_beta',
        'v_alpha',
        'j_alpha',
        'v_beta',
        'j_beta',
        'species',
        'source',
    ]

    pmhc_columns = [
        'id',
        'peptide',
        'mhc_pseudo',
        'mhc1',
        'mhc2',
        'mhc_type',
        'species',
        'source',
    ]

    overlap_columns = [col for col in tcr_columns if col in pmhc_columns]
    tcr_rename = {col: 'tcr_' + col for col in overlap_columns}
    pmhc_rename = {col: 'pmhc_' + col for col in overlap_columns}

    background = []

    for _, sequence in sequences.iterrows():
        out_group = sequences[
            (sequences['peptide'] != sequence['peptide']) & (sequences['collated_cdrs'] != sequence['collated_cdrs'])
        ]

        for _ in range(background_ratio):
            random_seq = out_group.sample(1, random_state=random_state).iloc[0]

            tcr = sequence[tcr_columns].rename(tcr_rename)
            pmhc = random_seq[pmhc_columns].rename(pmhc_rename)

            id_ = tcr['tcr_id'] + '_' + pmhc['pmhc_id']

            background.append(dict(id=id_, **tcr.to_dict(), **pmhc.to_dict()))

    return pd.DataFrame(background)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    if args.seed:
        logger.info('Setting random seed to %d', args.seed)
        rng = np.random.default_rng(args.seed)

    else:
        rng = np.random.default_rng()

    logger.info('Loading data from %s', args.input)
    sequences = pd.read_csv(args.input)

    logger.info('Processing column names')
    sequences.insert(0, 'id', sequences.index.map(str))
    for column_name in 'id', 'species', 'source':
        sequences[f'tcr_{column_name}'] = sequences[column_name]
        sequences[f'pmhc_{column_name}'] = sequences[column_name]

    if args.generate_background:
        logger.info(
            'Generating background by shuffling TCR and pMHC pairs at a ratio of 1:%d (positives:negatives)',
            args.background_proportion,
        )
        sequences['background'] = False

        background_data = generate_background(sequences, args.background_proportion, random_state=rng)
        background_data['background'] = True

        sequences = pd.concat([sequences, background_data]).reset_index(drop=True)

    logger.info('Sampling %d samples', args.sample_size)
    sample = sequences.sample(args.sample_size, random_state=rng)

    logger.info('Outputting to %s', args.output)
    sample.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
