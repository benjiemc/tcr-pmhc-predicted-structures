from unittest import TestCase

import numpy as np
import pandas as pd

from tcr_pmhc_predicted_structures.utils import assign_mhc_class, find_common_groups, get_cdr_sequences


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


class TestAssignMHCClass(TestCase):
    def test_class_i_human(self):
        self.assertEqual(assign_mhc_class('HLA-A*02:01'), 'MH1')

    def test_class_ii_human(self):
        self.assertEqual(assign_mhc_class('HLA-DQ*02:01'), 'MH2')

    def test_class_i_mouse(self):
        self.assertEqual(assign_mhc_class('H2-D'), 'MH1')

    def test_class_ii_mouse(self):
        self.assertEqual(assign_mhc_class('H2-IA'), 'MH2')


class TestGetCDRSequences(TestCase):
    def test_alpha_chain(self):
        cdr1_alpha, cdr2_alpha, cdr3_alpha = get_cdr_sequences(
            'MDAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASTWGRASTDTQY'
            'FGPGTRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRC'
            'QVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD',
        )

        self.assertEqual(cdr1_alpha, 'SGHNS')
        self.assertEqual(cdr2_alpha, 'FNNNVP')
        self.assertEqual(cdr3_alpha, 'ASTWGRASTDTQY')

    def test_beta_chain(self):
        cdr1_beta, cdr2_beta, cdr3_beta = get_cdr_sequences(
            'MAQTVTQSQPEMSVQEAETVTLSCTYDTSESDYYLFWYKQPPSRQMILVIRQEAYKQQNATENRFSVNFQKAAKSFSLKISDSQLGDAAMYFCASSGNTPLVFGKG'
            'TRLSVIPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS'
        )

        self.assertEqual(cdr1_beta, 'TSESDYY')
        self.assertEqual(cdr2_beta, 'QEAYKQQN')
        self.assertEqual(cdr3_beta, 'ASSGNTPLV')

    def test_invalid(self):
        cdr1, cdr2, cdr3 = get_cdr_sequences(
            'CLFTDFDSQTNVSQSKDSDVNSIIPEDTRLSVIPNSMDFKSNSAVAWSNKSDFACANAIQNPDPAVYQLRDSKSSYITDKCVLDMRFNDKSVTFFPSPESSMAQTV'
            'TQSQPEMLVIRQEAYKQQNATENRFSMYFCSVQEAETVTLSCTYDTSVNFQKAAKSFSLKISDSQLGDAAESDYYLFWYKQPPSRQMIASSGNTPLVFGKG'
        )

        self.assertIsNone(cdr1)
        self.assertIsNone(cdr2)
        self.assertIsNone(cdr3)
