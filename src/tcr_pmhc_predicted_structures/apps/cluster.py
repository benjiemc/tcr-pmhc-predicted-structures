import argparse
from warnings import simplefilter

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import ClusterWarning

parser = argparse.ArgumentParser(prog='ClusterSequences')

parser.add_argument('input', nargs='+')
parser.add_argument('--output', '-o', required=True)
parser.add_argument('--number-of-clusters', '-n', default=1000, type=int)
parser.add_argument('--normalize', action='store_true')


def main():
    args = parser.parse_args()

    data = [np.loadtxt(input_file) for input_file in args.input]

    composite = np.zeros(data[0].shape)

    for layer in data:
        if args.normalize:
            layer = (layer - layer.mean()) / layer.std()

        composite += layer

    simplefilter("ignore", ClusterWarning)

    clusters = AgglomerativeClustering(n_clusters=args.number_of_clusters,
                                       linkage='complete',
                                       affinity='precomputed').fit_predict(composite)

    np.savetxt(args.output, clusters.astype(int), fmt='%i')


if __name__ == '__main__':
    main()
