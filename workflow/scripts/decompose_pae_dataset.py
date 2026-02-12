#!/usr/bin/env python
"""Decompose PAE dataset using Principle Component Analysis."""

import logging

import h5py
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])

logger.info('Loading PAE dataset from %s', snakemake.input[0])
with h5py.File(snakemake.input[0], 'r') as fh:
    paes = fh['pae'][:]

    names = np.array([name.decode() for name in fh['id'][:]])
    seeds = fh['seed'][:]
    samples = fh['sample'][:]

    folds = fh['fold'][:]

logger.info('Decomposing dataset into %d component', snakemake.params.num_components)
pca = PCA(n_components=snakemake.params.num_components).fit_transform(paes)

pca_df = pd.DataFrame(
    {
        'name': names,
        'seed': seeds,
        'sample': samples,
        'pc1': pca[:, 0],
        'pc2': pca[:, 1],
    }
)

logger.info('Saving PCA components to %s', snakemake.output[0])
pca_df.to_csv(snakemake.output[0], index=False)
