from unittest import TestCase

import numpy as np
import pandas as pd

from tcr_pmhc_predicted_structures.utils import find_common_groups


class TestFindCommonGroups(TestCase):
    def test_none(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        with self.assertRaises(ValueError):
            find_common_groups(data, [])

    def test_1d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 2, 3, 3, 3, 4, 4, 5, 5]))

    def test_2d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_b': ['l', 'm', 'l', 'o', 'p', 'q', 'q', 'r', 's', 't'],
                'col_c': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a', 'col_b'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 1, 2, 2, 2, 2, 2, 3, 3]))

    def test_4d(self):
        data = pd.DataFrame(
            {
                'col_a': ['a', 'a', 'b', 'c', 'c', 'c', 'd', 'd', 'e', 'e'],
                'col_b': ['l', 'm', 'l', 'o', 'p', 'q', 'q', 'r', 's', 't'],
                'col_c': ['x', 'y', 'z', 'aa', 'bb', 'cc', 'dd', 'ee', 'ee', 'ee'],
                'col_d': ['aaa', 'bbb', 'ccc', 'ddd', 'eee', 'fff', 'ggg', 'ggg', 'ggg', 'ggg'],
                'col_e': [1, 1, 2, 3, 0, 1, 5, 6, 1, 10],
            }
        )

        groups = find_common_groups(data, ['col_a', 'col_b', 'col_c', 'col_d'])
        np.testing.assert_array_equal(groups, np.array([1, 1, 1, 2, 2, 2, 2, 2, 2, 2]))
