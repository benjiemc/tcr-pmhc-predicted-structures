"""Compute distances."""

import os
import sys

import numpy as np
import pandas as pd
from pyxdameraulevenshtein import damerau_levenshtein_distance


def levenshtein_distance(token1, token2):
    """Calculate Levenshtein distances.

    WARNING/TODO: Replaced by a Cython implementation.

    """
    distances = np.zeros((len(token1) + 1, len(token2) + 1))

    for t1 in range(len(token1) + 1):
        distances[t1][0] = t1

    for t2 in range(len(token2) + 1):
        distances[0][t2] = t2

    a = 0
    b = 0
    c = 0

    for t1 in range(1, len(token1) + 1):
        for t2 in range(1, len(token2) + 1):
            if token1[t1 - 1] == token2[t2 - 1]:
                distances[t1][t2] = distances[t1 - 1][t2 - 1]
            else:
                a = distances[t1][t2 - 1]
                b = distances[t1 - 1][t2]
                c = distances[t1 - 1][t2 - 1]

                if a <= b and a <= c:
                    distances[t1][t2] = a + 1
                elif b <= a and b <= c:
                    distances[t1][t2] = b + 1
                else:
                    distances[t1][t2] = c + 1

    return distances[len(token1)][len(token2)]


def create_distance_matrix(sequences):
    distance_matrix = np.zeros((len(sequences), len(sequences)))

    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i != j:
                distance_matrix[j][i] = damerau_levenshtein_distance(seq1, seq2)

    return distance_matrix


def main():
    input_data = sys.argv[1]
    name = sys.argv[2]
    output_dir = sys.argv[3]

    vdj_db_sequences = pd.read_csv(input_data)
    sequences = vdj_db_sequences[name].to_list()

    matrix = create_distance_matrix(sequences)

    np.savetxt(os.path.join(output_dir, name + '_distance_matrix.txt'), matrix)


if __name__ == '__main__':
    main()
