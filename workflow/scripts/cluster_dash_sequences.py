#!/usr/bin/env python
"""Cluster dash sequences."""

import logging

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.clustering import calculate_sequence_distances

logger = logging.getLogger()
setup_logger(logger, snakemake.config['log_level'])

logger.info('Loading data from %s', snakemake.input[0])
dash = pd.read_csv(snakemake.input[0])
dash.index.name = 'id'
dash = dash.reset_index()
dash['id'] = 'dash_' + dash['id'].astype(str)

logger.info('Dropping sequences with NaNs')
dash = dash.dropna()

if snakemake.wildcards.epitope != 'all':
    logger.info('Selecting TCRs binding to the %s epitope', snakemake.wildcards.epitope)
    dash = dash[dash['epitope'] == snakemake.wildcards.epitope].copy()

logger.info('Computing TCR levenshtein distances')
distance_matrix = calculate_sequence_distances(dash['cdr3_a_aa'].tolist(), dash['cdr3_b_aa'].tolist())
cdr3_distance_matrix = distance_matrix + distance_matrix.T

logger.info('Clustering CDR3 distances at a Levenshtein distance of %s', snakemake.params.distance_threshold)
clusters = AgglomerativeClustering(
    n_clusters=None,
    distance_threshold=snakemake.params.distance_threshold,
    compute_full_tree=True,
    linkage='complete',
    metric='precomputed',
).fit_predict(cdr3_distance_matrix)

logger.info('Filtering clusters smaller than %s', snakemake.params.min_cluster_size)
cluster_counts = {}
for cluster in clusters:
    if cluster not in cluster_counts:
        cluster_counts[cluster] = 1

    else:
        cluster_counts[cluster] += 1

clusters_filtered = np.array(
    [str(cluster) if cluster_counts[cluster] >= snakemake.params.min_cluster_size else 'None' for cluster in clusters]
)

logger.info('Writing clusters to %s', snakemake.output[0])
with open(snakemake.output[0], 'w') as fh:
    for id_, cluster in zip(dash['id'].tolist(), clusters_filtered, strict=True):
        fh.write(f'{id_}\t{cluster}\n')
