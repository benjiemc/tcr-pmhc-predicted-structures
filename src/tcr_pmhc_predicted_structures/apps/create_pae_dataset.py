"""Create an h5 dataset from PAE matrices."""

import argparse
import json
import logging
import os
import re
import sys

import h5py
import numpy as np
import pandas as pd

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    '--format',
    choices=['txt', 'json'],
    default='txt',
    help="format the PAE confidence scores are stored in, either 'txt' or 'json' (Default: 'txt')",
)
parser.add_argument(
    '--include-training-metadata',
    action='store_true',
    help='include label, fold, and source in final dataset',
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument(
    '--sequences', required=True, help='path to csv containing the sequences to make into a PAE dataset'
)
inputs.add_argument('--confidences', nargs='+', required=True, help='path to the PAE files formated as txt files')

outputs = parser.add_argument_group('Outputs')
outputs.add_argument('--output', '-o', required=True, help='path to output h5 dataset')

add_logging_arguments(parser)


def load_pae(path: str, format_: str) -> np.ndarray:
    match format_:
        case 'txt':
            pae = np.loadtxt(path, dtype=np.float32)

        case 'json':
            with open(path, 'r') as fh:
                confidences = json.load(fh)

            pae = np.array(confidences['pae'], dtype=np.float32)

    return pae


def main():
    args = parser.parse_args()
    setup_logger(logger, level=args.log_level)

    logger.info('Getting max size and max value')
    num_confidences = len(args.confidences)

    max_num_rows = 0
    max_num_cols = 0
    max_value = -np.inf
    for num, path in enumerate(args.confidences, 1):
        logger.debug('Loading confidence scores from %s - %d of %d', path, num, num_confidences)
        pae = load_pae(path, args.format)

        max_num_rows = max(max_num_rows, pae.shape[0])
        max_num_cols = max(max_num_cols, pae.shape[1])
        max_value = max(max_value, np.nanmax(pae))

    logger.debug('Max shape is (%d, %d), max value is %f', max_num_rows, max_num_cols, max_value)

    logger.info('Loading sequences from %s', args.sequences)
    sequences = pd.read_csv(args.sequences, index_col='id')
    sequences.index = sequences.index.astype('str')

    logger.info('Creating dataset at %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        str_dt = h5py.string_dtype(encoding='utf-8')
        ids = fh.create_dataset('id', shape=(num_confidences,), dtype=str_dt)
        seeds = fh.create_dataset('seed', shape=(num_confidences,), dtype=np.int32)
        samples = fh.create_dataset('sample', shape=(num_confidences,), dtype=np.int32)

        if args.include_training_metadata:
            sources = fh.create_dataset('source', shape=(num_confidences,), dtype=np.int8)
            folds = fh.create_dataset('fold', shape=(num_confidences,), dtype=np.int32)
            labels = fh.create_dataset('label', shape=(num_confidences,), dtype=np.int8)

        paes = fh.create_dataset(
            'pae',
            shape=(num_confidences, max_num_rows * max_num_cols),
            dtype=np.float32,
            chunks=(1, max_num_rows * max_num_cols),
            compression='gzip',
            compression_opts=4,
            fillvalue=0.0,
        )

        for i, path in enumerate(args.confidences):
            logger.debug('Adding %d of %d', i + 1, num_confidences)

            entry_name = re.search(
                r'^(?P<name>[A-Za-z0-9_]+)_(?:seed-)?(?P<seed>\d+)_(?:sample-)?(?P<sample>\d+)',
                os.path.basename(path),
            ).groupdict()

            logger.debug('Collecting info from %s in sequences', entry_name['name'])
            info = sequences.loc[entry_name['name']]

            logger.debug('Loading confidence scores from %s', path)
            pae = load_pae(path, args.format)

            logger.debug('Padding PAE')
            pae = np.pad(
                pae,
                ((0, max_num_rows - pae.shape[0]), (0, max_num_cols - pae.shape[1])),
                mode='constant',
                constant_values=np.nan,
            )

            logger.debug('Inverting and normalising the PAE')
            pae = (max_value - pae) / max_value

            logger.debug('Filling NaNs as zeros')
            pae = np.nan_to_num(pae, nan=0.0)

            logger.debug('Flattening PAE matrices')
            pae = pae.flatten()

            logger.debug('Writing line')
            ids[i] = entry_name['name']
            seeds[i] = int(entry_name['seed'])
            samples[i] = int(entry_name['sample'])

            if args.include_training_metadata:
                sources[i] = int(info['source'] == 'Sequences')
                folds[i] = info['fold']
                labels[i] = info['label']

            paes[i, :] = pae


if __name__ == '__main__':
    main()
