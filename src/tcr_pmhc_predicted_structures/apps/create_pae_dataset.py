"""Create an h5 dataset from PAE matrices."""

import argparse
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

inputs = parser.add_argument_group('Inputs')
inputs.add_argument(
    '--sequences', required=True, help='path to csv containing the sequences to make into a PAE dataset'
)
inputs.add_argument('--confidences', nargs='+', required=True, help='path to the PAE files formated as txt files')

outputs = parser.add_argument_group('Outputs')
outputs.add_argument('--output', '-o', required=True, help='path to output h5 dataset')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, level=args.log_level)

    logger.info('Loading sequences from %s', args.sequences)
    sequences = pd.read_csv(args.sequences, index_col='id')
    sequences.index = sequences.index.astype('str')

    logger.info('Loading PAEs')
    ids = []
    seeds = []
    samples = []
    sources = []
    paes = []
    folds = []
    labels = []
    for path in args.confidences:
        logger.debug('Loading PAE from %s', path)
        basename = os.path.basename(path)

        pae = np.loadtxt(path)
        pae = pae.flatten()

        entry_name = re.search(r'^(?P<name>[A-Za-z0-9_]+)_(?P<seed>\d+)_(?P<sample>\d+)', basename).groupdict()

        logger.debug('Collecting info from %s in sequences', entry_name['name'])
        info = sequences.loc[entry_name['name']]

        ids.append(entry_name['name'])
        seeds.append(int(entry_name['seed']))
        samples.append(int(entry_name['sample']))

        sources.append(int(info['source'] == 'Sequences'))
        folds.append(info['fold'])
        labels.append(info['label'])

        paes.append(pae)

    seeds = np.array(seeds)
    samples = np.array(samples)
    folds = np.array(folds)
    labels = np.array(labels)

    paes = np.stack(paes)

    logger.info('Inverting and normalising the PAEs')
    max_value = np.nanmax(paes)
    paes = (max_value - paes) / max_value

    logger.info('Filling NaNs as zeros')
    paes = np.nan_to_num(paes)

    logger.info('Writing dataset to %s', args.output)
    with h5py.File(args.output, 'w') as fh:
        fh['id'] = ids
        fh['seed'] = seeds
        fh['sample'] = samples
        fh['pae'] = paes
        fh['source'] = sources
        fh['fold'] = folds
        fh['label'] = labels


if __name__ == '__main__':
    main()
