"""Compute averages."""

import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('--clusters', help='path to cluster file')
parser.add_argument('--distance-matrices', nargs='+', help='paths to distance matrices')
parser.add_argument('--output', help='path to output csv')


if __name__ == '__main__':
    args = parser.parse_args()

    print('Loading data...')
    clusters = np.loadtxt(args.clusters).astype(int)
    num_points = len(clusters)

    combo = np.zeros((num_points, num_points))
    matrices = []

    for path in args.distance_matrices:
        name = path.split('/')[-1].split('.')[0]
        distance_matrix = np.loadtxt(path).astype(int)
        matrices.append((name, distance_matrix))

        combo += distance_matrix

    matrices.append(('combination', combo))

    cluster_names = []
    distance_names = []
    in_avgs = []
    out_avgs = []

    for cluster in np.unique(clusters):
        print('Working on cluster', cluster)

        in_cluster_idx = np.where(cluster == clusters, 1, 0)
        out_cluster_idx = np.where(cluster != clusters, 1, 0)

        in_cluster_mask = np.outer(in_cluster_idx, in_cluster_idx)
        out_cluster_mask = np.outer(out_cluster_idx, out_cluster_idx)

        for name, distances in matrices:
            in_cluster_dist = np.where(in_cluster_mask, distances, np.nan)
            out_cluster_dist = np.where(out_cluster_mask, distances, np.nan)

            in_avgs.append(np.nanmean(in_cluster_dist))
            out_avgs.append(np.nanmean(out_cluster_dist))

            cluster_names.append(cluster)
            distance_names.append(name)

    print('Writing data to file')
    cluster_distances = pd.DataFrame(
        {
            'cluster': cluster_names,
            'distance_name': distance_names,
            'in_avg_distance': in_avgs,
            'out_avg_distance': out_avgs,
        }
    )

    cluster_distances.to_csv(args.output, index=False)
