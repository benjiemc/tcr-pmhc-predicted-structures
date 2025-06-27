"""Create contact map between a TCR and pMHC."""

import argparse
import logging
import re
import sys

import numpy as np
import pandas as pd
import stcrpy

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, assign_cdr_number
from tcr_pmhc_predicted_structures.utils import stcrpy_to_pandas

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument('structure', help='path to the input structure (mmCIF format)')

output = parser.add_argument_group('Output')
output.add_argument('--output', '-o', required=True, help='path to output CSV')

parser.add_argument(
    '--cutoff-distance',
    default=5.0,
    type=float,
    help='maximum distance to consider a contact (default: 5.0 Å)',
)
parser.add_argument(
    '--even-offset-side',
    choices=['left', 'right'],
    default='left',
    help="Offset side when relative_pos_centre is selected (Default: 'left')",
)

add_logging_arguments(parser)


HEAVY_ATOMS: set[str] = {'C', 'N', 'O', 'S'}
'''Heavy atoms that make up amino acids.'''

CHAIN_ANNOTATIONS = {
    'A': 'mhc_chain1',
    'B': 'mhc_chain2',
    'C': 'antigen_chain',
    'D': 'Achain',
    'E': 'Bchain',
}
'''Annotations for standardised chain names.'''


def enumerate_chain(chain_df: pd.DataFrame) -> pd.Series:
    """Enumerate a chain by residue to get a residue number."""
    chain_df = chain_df.copy()
    chain_df['residue_insert_code'] = chain_df['residue_insert_code'].fillna('')

    return chain_df.groupby(['residue_seq_id', 'residue_insert_code']).ngroup() + 1


def number_relative_to_centre(numbering: np.array, even_offset='left') -> np.array:
    """Numbers a sequential array relative to its centre.

    Examples:
        # Odd length example
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6, 7]))
        array([-3, -2, -1, 0, 1, 2, 3])

        # Even length example - left offset
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6]), even_offset='left')
        array([-2, -1, 0, 1, 2, 3])

        # Even length example - right offset
        >>> number_relative_to_centre(np.array([1, 2, 3, 4, 5, 6]), even_offset='left')
        array([-3, -2, -1, 0, 1, 2])

    Args:
        numbering: sequential input array
        even_offset: offset to give even length input numbering, can be either 'left' or 'right' (Default: 'left')

    Returns:
        array of the same length but with numbering relative to centre.

    Raises:
        ValueError: if even_offset is not 'left' or 'right'

    """
    length = np.max(numbering)

    if length % 2 == 0:
        match even_offset:
            case 'left':
                return numbering - length // 2

            case 'right':
                return numbering - (length // 2 + 1)

            case _:
                msg = f"Invalid even_offset={even_offset}. Must be 'left' or 'right'."
                raise ValueError(msg)

    return numbering - np.ceil(length / 2)


def collect_contacts(
    structure_df: pd.DataFrame,
    mhc_type: str,
    cutoff_distance: float,
    even_offset_side: str = 'left',
) -> pd.DataFrame:
    """Collect the contacting residues from the structure."""
    logger.debug('Annotating structure')
    structure_df['resi'] = structure_df['residue_seq_id'].apply(str) + structure_df['residue_insert_code'].fillna('')
    structure_df['chain_type'] = structure_df['chain_id'].map(CHAIN_ANNOTATIONS)
    structure_df['cdr'] = None
    structure_df.loc[
        (structure_df['chain_type'] == 'Achain') | (structure_df['chain_type'] == 'Bchain'),
        'cdr',
    ] = structure_df[(structure_df['chain_type'] == 'Achain') | (structure_df['chain_type'] == 'Bchain')][
        'residue_seq_id'
    ].map(assign_cdr_number)

    logger.debug('Splitting TCR and pMHC')
    tcr_cdrs_df = structure_df[structure_df['cdr'].notna() & structure_df['element'].isin(HEAVY_ATOMS)].copy()

    tcr_cdrs_df['cdr_position'] = (
        tcr_cdrs_df.groupby(['chain_type', 'cdr'])
        .apply(
            enumerate_chain,
            include_groups=False,
        )
        .droplevel(['chain_type', 'cdr'])
    )

    tcr_cdrs_df['relative_pos_centre'] = tcr_cdrs_df.groupby(
        ['chain_type', 'cdr'],
    )['cdr_position'].transform(number_relative_to_centre, even_offset=even_offset_side)

    match mhc_type:
        case 'MH1':
            mhc_abd_df = structure_df[
                (structure_df['chain_type'] == 'mhc_chain1')
                & structure_df['residue_seq_id'].isin(IMGT_MH1_ABD)
                & structure_df['element'].isin(HEAVY_ATOMS)
            ].copy()

        case 'MH2':
            mhc_abd_df = structure_df[
                ((structure_df['chain_type'] == 'mhc_chain1') | (structure_df['chain_type'] == 'mhc_chain2'))
                & structure_df['residue_seq_id'].isin(IMGT_MH2_ABD)
                & structure_df['element'].isin(HEAVY_ATOMS)
            ].copy()

        case _:
            msg = f'Invalid MHC Type: {mhc_type}'
            raise ValueError(msg)

    peptide_df = structure_df[
        (structure_df['chain_type'] == 'antigen_chain') & structure_df['element'].isin(HEAVY_ATOMS)
    ].copy()

    peptide_df['peptide_position'] = enumerate_chain(peptide_df)
    peptide_df['relative_pos_centre'] = number_relative_to_centre(
        peptide_df['peptide_position'],
        even_offset=even_offset_side,
    )
    peptide_df['resi'] = None

    pmhc_df = pd.concat([peptide_df, mhc_abd_df])

    logger.debug('Calculating distances')
    interface = tcr_cdrs_df.merge(pmhc_df, how='cross', suffixes=('_tcr', '_pmhc'))
    interface['distance'] = np.sqrt(
        np.square(interface['pos_x_tcr'] - interface['pos_x_pmhc'])
        + np.square(interface['pos_y_tcr'] - interface['pos_y_pmhc'])
        + np.square(interface['pos_z_tcr'] - interface['pos_z_pmhc'])
    )

    logger.debug('Getting atoms within %.2f Å of each other', cutoff_distance)
    interface_contacts = interface[interface['distance'] < cutoff_distance].copy()
    interface_contacts['mhc_type'] = mhc_type

    return interface_contacts


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading structure from %s', args.structure)
    tcr = stcrpy.load_TCR(args.structure)

    tcr_df = stcrpy_to_pandas(tcr)
    mhc_df = stcrpy_to_pandas(tcr.get_MHC()[0])
    antigen_df = stcrpy_to_pandas(tcr.antigen)

    structure_df = pd.concat([tcr_df, mhc_df, antigen_df]).reset_index(drop=True)

    logger.info('Collecting contacts')
    contacts = collect_contacts(
        structure_df,
        re.search(r'MH\d', str(type(tcr.get_MHC()[0])))[0],
        args.cutoff_distance,
        even_offset_side=args.even_offset_side,
    )

    logger.info('Normalising outputs')
    contacts = contacts.drop_duplicates(
        [
            'chain_type_tcr',
            'residue_seq_id_tcr',
            'residue_insert_code_tcr',
            'chain_type_pmhc',
            'residue_seq_id_pmhc',
            'residue_insert_code_pmhc',
        ]
    )
    contacts['cdr_type'] = (
        'CDR'
        + contacts['cdr_tcr'].apply(int).apply(str)
        + contacts['chain_type_tcr'].map({'Achain': 'alpha', 'Bchain': 'beta'})
    )

    output_columns = [
        'cdr_type',
        'chain_id_tcr',
        'chain_type_tcr',
        'residue_seq_id_tcr',
        'residue_insert_code_tcr',
        'resi_tcr',
        'residue_name_tcr',
        'cdr_tcr',
        'cdr_position',
        'relative_pos_centre_tcr',
        'chain_id_pmhc',
        'residue_seq_id_pmhc',
        'residue_insert_code_pmhc',
        'resi_pmhc',
        'residue_name_pmhc',
        'chain_type_pmhc',
        'peptide_position',
        'relative_pos_centre_pmhc',
        'mhc_type',
    ]
    contacts = contacts[output_columns]
    contacts = contacts.sort_values(output_columns).reset_index(drop=True)

    logger.info('Outputting results to %s', args.output)
    contacts.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
