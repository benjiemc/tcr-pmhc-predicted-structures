from unittest import TestCase

import numpy as np

from tcr_pmhc_predicted_structures.clustering import calculate_sequence_distances


class TestCalculateSequenceDistances(TestCase):
    def test(self):
        alpha_seqs = ['CAB', 'CAAB', 'TASK']
        beta_seqs = ['CAT', 'BAT', 'YES']
        distance_matrix = calculate_sequence_distances(alpha_seqs, beta_seqs)

        np.testing.assert_array_equal(distance_matrix, np.array([[0, 1, 3], [1, 0, 3], [3, 3, 0]]))
