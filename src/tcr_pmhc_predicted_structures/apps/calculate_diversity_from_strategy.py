"""Calculate the diversity of a sample based on the sampling strategy."""

import argparse
import logging
import re

import numpy as np
import pandas as pd

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.diversity_metrics import pielou_entropy, shannon_entropy, simpson_index

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser()

parser.add_argument('--num-samples', '-n', type=int, default=50, help='Number of samples to take')
parser.add_argument('--sequences', help='Path to sequence file')
parser.add_argument('--clusters', help='Cluster annotations matching the number of the sequences')
parser.add_argument(
    '--baseline', action='store_true', help='Use uniform sampling as a baseline. clusters are note required.'
)
parser.add_argument('--num-sequences', type=int, help='Number of sequences to take in a sample')
parser.add_argument('--output', '-o', help='Output path')


def convert_allele_to_slug(code: str) -> str:
    """Convert allele to slug, eg. HLA-A*02:01:59 -> hla_a_02_01_59."""
    slug = code.lower()
    slug = re.sub(r'[*:-]', '_', slug)

    return slug


def main():
    setup_logger(logger, level='debug')
    args = parser.parse_args()

    sequences = pd.read_csv(args.sequences)

    if args.clusters:
        sequences['cluster'] = np.loadtxt(args.clusters).astype(int).astype(str)

    mhc_slugs = sequences['mhc_a_allele'].apply(convert_allele_to_slug)
    mhc_slugs_cropped = mhc_slugs.map(lambda slug: '_'.join(slug.split('_')[:4]))
    sequences['mhc_slugs'] = mhc_slugs_cropped

    categories = []
    shannons = []
    pielous = []
    simpsons = []

    for _ in range(args.num_samples):
        if args.baseline:
            logger.info('Generating sample from random sampling')
            sample = sequences.sample(args.num_sequences)

        else:
            logger.info('Generating sample from random sampling of clusters')
            sample = sequences.groupby('cluster').sample(1)

        sample_size = len(sample)

        for category in (
            'cdr3_alpha',
            'cdr3_beta',
            'epitope',
            'mhc_slugs',
            'mhc_pseudo_sequence',
            'tcr_pseudo_sequence',
        ):
            probabilities = sample.value_counts(category).to_numpy() / sample_size

            shannon = shannon_entropy(probabilities)
            pielou = pielou_entropy(probabilities)
            simpson = simpson_index(probabilities)

            categories.append(category)
            shannons.append(shannon)
            pielous.append(pielou)
            simpsons.append(simpson)

    logger.info('Exporting to csv file')
    pd.DataFrame(
        {
            'category': categories,
            'shannon_entropy': shannons,
            'pielou_entropy': pielous,
            'simpson': simpsons,
        }
    ).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
