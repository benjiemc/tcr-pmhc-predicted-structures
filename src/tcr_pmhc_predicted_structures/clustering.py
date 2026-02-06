"""Functions for clustering TCRs and pMHCs."""

import numpy as np
from pyxdameraulevenshtein import damerau_levenshtein_distance


def calculate_sequence_distances(
    cdr_3_alpha_sequences: list[str], cdr_3_beta_sequences: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the levenshtein distances between CDR3 loops.

    Args:
        cdr_3_alpha_sequences: list of CDR3 alpha sequences, must be the same length as the cdr_3_beta_sequences
        cdr_3_beta_sequences: list of CDR3 beta sequences, must be the same length as the cdr_3_alpha_sequences

    Returns:
        an asymetric distance matrix with the alpha sequence distances in the upper diagonal and the beta sequences in
        the lower diagonal

    """
    num_tcrs = len(cdr_3_alpha_sequences)

    cdr3_alpha_distance_matrix = np.zeros((num_tcrs, num_tcrs), dtype=int)
    for j in range(num_tcrs - 1):
        for i in range(j + 1, num_tcrs):
            cdr3_alpha_distance_matrix[i, j] = damerau_levenshtein_distance(
                cdr_3_alpha_sequences[i], cdr_3_alpha_sequences[j]
            )

    cdr3_beta_distance_matrix = np.zeros((num_tcrs, num_tcrs), dtype=int)
    for i in range(num_tcrs - 1):
        for j in range(i + 1, num_tcrs):
            cdr3_beta_distance_matrix[i, j] = damerau_levenshtein_distance(
                cdr_3_beta_sequences[i], cdr_3_beta_sequences[j]
            )

    return np.max([cdr3_alpha_distance_matrix, cdr3_beta_distance_matrix], axis=0)
