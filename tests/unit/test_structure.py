import os
from unittest import TestCase

import pandas as pd
from Bio.PDB import PDBParser

from tcr_pmhc_predicted_structures.structure import bio_to_pandas, extract_chains, get_sequence

TEST_DATA = 'tests/data'


class TestGetSequence(TestCase):
    def setUp(self):
        parser = PDBParser()
        self.structure = parser.get_structure('test', os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

    def test(self):
        sequence = get_sequence(self.structure, 'D')
        self.assertEqual(
            sequence,
            (
                'KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSD'
                'SATYLCAVTTDSWGKLQFGAGTQVVVTP'
            ),
        )

    def test_res_selection(self):
        sequence = get_sequence(self.structure, 'D', residue_range=set(range(105, 117 + 1)))
        self.assertEqual(sequence, ('AVTTDSWGKLQ'))


class TestExtractChains(TestCase):
    def setUp(self):
        parser = PDBParser()
        self.structure = parser.get_structure('test', os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

    def test_e(self):
        test_struct = extract_chains(self.structure, ['E'])
        test_chains = {chain.id for chain in test_struct[0]}

        self.assertEqual(test_chains, {'E'})

    def test_de(self):
        test_struct = extract_chains(self.structure, ['D', 'E'])
        test_chains = {chain.id for chain in test_struct[0]}

        self.assertEqual(test_chains, {'D', 'E'})

    def test_deca(self):
        test_struct = extract_chains(self.structure, ['D', 'E', 'C', 'A'])
        test_chains = {chain.id for chain in test_struct[0]}

        self.assertEqual(test_chains, {'D', 'E', 'C', 'A'})


class TestBioToPandas(TestCase):
    def test(self):
        parser = PDBParser()
        structure = parser.get_structure('test', os.path.join(TEST_DATA, '1qse_DECAB.pdb'))

        structure_df = bio_to_pandas(structure)

        columns = [
            'record_type',
            'atom_number',
            'atom_name',
            'alt_loc',
            'residue_name',
            'chain_id',
            'residue_seq_id',
            'residue_insert_code',
            'pos_x',
            'pos_y',
            'pos_z',
            'occupancy',
            'b_factor',
            'element',
            'charge',
        ]

        pd.testing.assert_frame_equal(
            structure_df.head(),
            pd.DataFrame(
                [
                    ['ATOM', 1, 'N', None, 'LEU', 'C', 1, None, 80.148003, -9.695, 37.671001, 1.0, 5.00, 'N', None],
                    ['ATOM', 2, 'CA', None, 'LEU', 'C', 1, None, 80.480003, -10.378, 38.950001, 1.0, 17.98, 'C', None],
                    ['ATOM', 3, 'C', None, 'LEU', 'C', 1, None, 80.158997, -9.536, 40.179001, 1.0, 19.34, 'C', None],
                    ['ATOM', 4, 'O', None, 'LEU', 'C', 1, None, 80.869003, -8.570, 40.467999, 1.0, 20.10, 'O', None],
                    ['ATOM', 5, 'CB', None, 'LEU', 'C', 1, None, 81.962997, -10.682, 39.001999, 1.0, 22.45, 'C', None],
                ],
                columns=columns,
            ),
            check_dtype=False,
        )

        pd.testing.assert_frame_equal(
            structure_df.tail(),
            pd.DataFrame(
                [
                    [
                        'ATOM',
                        3284,
                        'CD',
                        None,
                        'ARG',
                        'A',
                        1091,
                        None,
                        74.278999,
                        -6.204,
                        17.808001,
                        1.0,
                        90.67,
                        'C',
                        None,
                    ],
                    [
                        'ATOM',
                        3285,
                        'NE',
                        None,
                        'ARG',
                        'A',
                        1091,
                        None,
                        73.511002,
                        -6.038,
                        19.045000,
                        1.0,
                        89.83,
                        'N',
                        None,
                    ],
                    [
                        'ATOM',
                        3286,
                        'CZ',
                        None,
                        'ARG',
                        'A',
                        1091,
                        None,
                        72.431000,
                        -6.761,
                        19.346001,
                        1.0,
                        87.80,
                        'C',
                        None,
                    ],
                    [
                        'ATOM',
                        3287,
                        'NH1',
                        None,
                        'ARG',
                        'A',
                        1091,
                        None,
                        71.989998,
                        -7.701,
                        18.511000,
                        1.0,
                        94.21,
                        'N',
                        None,
                    ],
                    [
                        'ATOM',
                        3288,
                        'NH2',
                        None,
                        'ARG',
                        'A',
                        1091,
                        None,
                        71.759003,
                        -6.550,
                        20.468000,
                        1.0,
                        91.82,
                        'N',
                        None,
                    ],
                ],
                columns=columns,
                index=[3283, 3284, 3285, 3286, 3287],
            ),
            check_dtype=False,
        )
