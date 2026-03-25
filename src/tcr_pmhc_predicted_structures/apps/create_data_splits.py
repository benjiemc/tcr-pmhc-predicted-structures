"""Create splits in the data and shuffle entries to create negative data.

Required columns:

    - cdr1_alpha
    - cdr2_alpha
    - cdr3_alpha
    - cdr1_beta
    - cdr2_beta
    - cdr3_beta
    - peptide
    - mhc_pseudo

"""

import argparse
import logging
import sys

import h5py
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.utils import create_even_folds, find_common_groups

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input_data', help='path to the input csv file')
parser.add_argument('--output', '-o', required=True, help='path to output csv file with the annotated sequences')
parser.add_argument('--seed', default=None, type=int, help='random seed for data splitting')

data_group = parser.add_argument_group('Data')
data_group.add_argument(
    '--negative-proportion',
    default=5,
    type=int,
    help=(
        'ratio of positives to negatives by randomly sampling the background for the negatives (Default: 5, meaning 5 '
        'times more negatives than positives)'
    ),
)

data_group.add_argument(
    '--num-folds',
    type=int,
    default=5,
    help='Number of folds to create in the dataset (Default: 5)',
)
data_group.add_argument(
    '--split-type',
    choices=['random', 'tcr', 'peptide', 'pMHC', 'levenshtein'],
    default='pMHC',
    help=(
        "Method to partition data between folds (Default: 'pMHC'). 'random' means to randomly shuffle data between"
        " folds, 'tcr' means no TCRs are shared across folds, 'peptide' means to ensure no peptides are shared across "
        "folds, 'pMHC' means that no peptides or MHCs (based on pseudo sequence) are shared across folds, and "
        "'levenshtein' means to use a levenshtein distance to separate data points between folds (more parameters "
        "below)."
    ),
)
data_group.add_argument(
    '--split-distance',
    default=6,
    type=int,
    help=(
        "If using a 'levenshtein' split, the mimimum distance between data points in two different folds (Default: 6)."
    ),
)
data_group.add_argument(
    '--split-entities',
    nargs='+',
    choices=[
        'cdr1_alpha',
        'cdr2_alpha',
        'cdr3_alpha',
        'cdr1_beta',
        'cdr2_beta',
        'cdr3_beta',
        'peptide',
        'mhc_pseudo',
    ],
    help="If using a 'levenshtein' split, the entities to consider as part of the split.",
)
data_group.add_argument(
    '--distances',
    help=(
        "If using a 'levenshtein' split, the path to an hdf5 file containing the distances between each type of entity "
        "considered."
    ),
)

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    if args.seed:
        logger.info('Setting seed to %d', args.seed)
        rng = np.random.default_rng(args.seed)

    else:
        rng = np.random.default_rng()

    logger.info('Loading data')
    sequence_data = pd.read_csv(args.input_data)
    sequence_data['label'] = 1

    logger.info('Generating negative data by random sampling')
    sequence_data['collated_cdrs'] = sequence_data.filter(
        regex=r'^cdr[1-3]_(alpha|beta)$',
    ).apply('-'.join, axis='columns')

    peptides = sequence_data['peptide'].unique()

    tcr_columns = sequence_data.filter(regex='cdr|(v|j)_(alpha|beta)|species').columns
    pmhc_columns = sequence_data.filter(regex='peptide|mhc|species').columns

    negative_data = []

    for peptide in peptides:
        in_group = sequence_data[sequence_data['peptide'] == peptide]
        out_group = sequence_data[
            (sequence_data['peptide'] != peptide) & (~sequence_data['collated_cdrs'].isin(in_group['collated_cdrs']))
        ]

        negatives = (
            out_group.sample(args.negative_proportion * len(in_group), random_state=rng)
            if len(out_group) >= (args.negative_proportion * len(in_group))
            else out_group
        )[tcr_columns]
        negatives = negatives.rename({'species': 'tcr_species'}, axis='columns')

        negatives = pd.concat(
            [
                negatives,
                pd.DataFrame(
                    [in_group[pmhc_columns].iloc[0].to_numpy()] * len(negatives),
                    columns=in_group[pmhc_columns].iloc[0].index,
                    index=negatives.index,
                ),
            ],
            axis=1,
        )
        negatives['species'] = negatives[['tcr_species', 'species']].apply(
            lambda row: '/'.join({row.tcr_species, row.species}), axis=1
        )
        negatives = negatives.drop('tcr_species', axis='columns')

        negatives['source'] = 'shuffling'
        negatives['label'] = 0

        negative_data.append(negatives)

    sequence_data = pd.concat([sequence_data, *negative_data]).reset_index(drop=True)

    logger.info('Creating %d cross-validation folds', args.num_folds)
    match args.split_type:
        case 'random':
            logger.debug('Randomly assigning cross-validation folds')

            if len(sequence_data) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d folds will be created',
                    args.num_folds,
                    len(sequence_data),
                )

            folds = create_even_folds(
                [(idx, 1) for idx in sequence_data.index],
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data.index.map({idx: i for i, fold in enumerate(folds, 1) for idx in fold})

        case 'tcr':
            logger.debug('Splitting TCRs across cross-validation folds')
            tcr_counts = sequence_data['collated_cdrs'].value_counts()
            tcr_counts = tcr_counts.sort_index().sort_values(ascending=False)

            if len(tcr_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(tcr_counts),
                )

            folds = create_even_folds(
                list(zip(tcr_counts.index.tolist(), tcr_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data['collated_cdrs'].map(
                {tcr_sequence: i for i, fold in enumerate(folds, 1) for tcr_sequence in fold}
            )

        case 'peptide':
            logger.debug('Splitting peptides across cross-validation folds')
            peptide_counts = sequence_data['peptide'].value_counts()
            peptide_counts = peptide_counts.sort_index().sort_values(ascending=False)

            if len(peptide_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(peptide_counts),
                )

            folds = create_even_folds(
                list(zip(peptide_counts.index.tolist(), peptide_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = sequence_data['peptide'].map(
                {peptide_sequence: i for i, fold in enumerate(folds, 1) for peptide_sequence in fold}
            )

        case 'pMHC':
            logger.debug('Splitting peptides across cross-validation folds')
            pmhc_groups = pd.Series(find_common_groups(sequence_data, ['peptide', 'mhc_pseudo']))
            pmhc_group_counts = pmhc_groups.value_counts()

            if len(pmhc_group_counts) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(pmhc_group_counts),
                )

            folds = create_even_folds(
                list(zip(pmhc_group_counts.index.tolist(), pmhc_group_counts.tolist(), strict=True)),
                num_folds=args.num_folds,
                seed=args.seed,
            )
            sequence_data['fold'] = pmhc_groups.map(
                {group_id: i for i, fold in enumerate(folds, 1) for group_id in fold}
            )

        case 'levenshtein':
            logger.debug(
                'Splitting data points across cross-validation folds based on a levenshtein distance of %d',
                args.split_distance,
            )

            with h5py.File(args.distances) as fh:
                distances = np.zeros((len(sequence_data), len(sequence_data)), dtype=int)

                for entity in args.split_entities:
                    names_map = {name.decode('utf-8'): idx for idx, name in enumerate(fh[entity]['names'][:])}
                    indices = sequence_data[entity].map(names_map).to_numpy()

                    distances += fh[entity]['distance_matrix'][:][np.ix_(indices, indices)]

            clusters = AgglomerativeClustering(
                n_clusters=None,
                distance_threshold=args.split_distance,
                metric='precomputed',
                linkage='single',
            ).fit_predict(distances)

            cluster_sizes = {}
            for cluster in clusters:
                if cluster in cluster_sizes:
                    cluster_sizes[cluster] += 1

                else:
                    cluster_sizes[cluster] = 1

            cluster_sizes = sorted(cluster_sizes.items(), key=lambda cluster_size: cluster_size[1], reverse=True)

            if len(cluster_sizes) < args.num_folds:
                logger.warning(
                    'Insufficient data to create %d folds. Only %d fold(s) will be created',
                    args.num_folds,
                    len(cluster_sizes),
                )

            folds = create_even_folds(cluster_sizes, num_folds=args.num_folds, seed=args.seed)
            fold_map = {cluster: i for i, fold in enumerate(folds, 1) for cluster in fold}
            sequence_data['fold'] = [fold_map[cluster] for cluster in clusters]

    sequence_data = sequence_data.drop('collated_cdrs', axis='columns')
    sequence_data.index.name = 'index'

    logger.info('Outputting data to %s', args.output)
    sequence_data.to_csv(args.output)


if __name__ == '__main__':
    main()
