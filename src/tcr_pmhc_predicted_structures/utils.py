"""Utility functions for processing TCR:pMHC data."""

import logging
import re
import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
from Bio.PDB import Entity
from Stitchr import stitchr as st
from Stitchr import stitchrfunctions as fxn

logger = logging.getLogger(__name__)


def mhc_code_to_slug(code: str) -> str:
    """Convert MHC allele codes into slugs.

    >>> mhc_code_to_slug('HLA-A*02:01:59')
    'hla_a_02_01_59'

    """
    slug = code.lower()
    slug = re.sub(r'[*:-]', '_', slug)

    return slug


def mhc_slug_to_code(slug: str) -> str:
    """Convert mhc slugs into the allele codes.

    >>> mhc_slug_to_code('hla_a_02_01')
    'HLA-A*02:01'

    >>> mhc_slug_to_code('h2_kb')
    'H2-Kb'

    TODO: Make this better for other mouse alleles
    """
    species = 'human' if slug.startswith('hla') else 'mouse' if slug.startswith('h2') else None

    if species == 'human':
        code = slug.upper()
        code = code.split('_')
        code = code[0] + '-' + code[1] + '*' + ':'.join(code[2:])

    elif species == 'mouse':
        code = slug.title()
        code = code.replace('_', '-')

    else:
        logger.error('Species not found, outputting slug')
        code = slug

    return code


def assign_species(allele: str) -> str | None:
    """Assign species based on MHC allele code."""
    if re.search('^HLA-', allele):
        return 'Human'

    if re.search('^H2-?', allele):
        return 'Mouse'

    if re.search('^Gaga', allele):
        return 'Chicken'

    return None


@lru_cache
def stitch_sequence(v_gene: str, j_gene: str, cdr3: str, species: str) -> str | None:
    """Create full length TCR sequence from V gene, J, gene, CDR3 and species information."""

    def create_input_args(args: dict, gene_types: list) -> tuple:
        input_args, chain = fxn.sort_input(args)
        codons = fxn.get_optimal_codons(input_args['codon_usage_path'], input_args['species'])
        j_res, low_conf_js = fxn.get_j_motifs(input_args['species'])
        c_res = fxn.get_c_motifs(input_args['species'])
        tcr_dat, tcr_functionality, partial = fxn.get_ref_data(chain, gene_types, input_args['species'])

        if input_args['extra_genes']:
            tcr_dat, tcr_functionality = fxn.get_additional_genes(tcr_dat, tcr_functionality)
            input_args['skip_c_checks'] = True

        if input_args['preferred_alleles_path']:
            preferred_alleles = fxn.get_preferred_alleles(
                input_args['preferred_alleles_path'],
                gene_types,
                tcr_dat,
                partial,
                chain,
            )

        else:
            preferred_alleles = {}

        return input_args, tcr_dat, tcr_functionality, partial, codons, preferred_alleles, c_res, j_res, low_conf_js

    gene_types = list(fxn.regions.values())

    start = 'C' if not cdr3.startswith('C') else ''
    end = 'F' if not cdr3.startswith('C') else ''

    cdr3 = start + cdr3 + end

    args = {
        'v': v_gene,
        'j': j_gene,
        'cdr3': cdr3,
        'species': species.upper(),
        'c': '',
        'l': '',
        'aa': '',
        'name': '',
        'seamless': False,
        '5_prime_seq': '',
        '3_prime_seq': '',
        'extra_genes': False,
        'mode': 'AA',
        'preferred_alleles_path': '',
        'codon_usage_path': '',
        'j_warning_threshold': 3,
        'skip_c_checks': False,
        'skip_n_checks': False,
        'suppress_warnings': False,
        'no_leader': True,
    }

    (
        input_args,
        tcr_dat,
        tcr_functionality,
        partial,
        codons,
        preferred_alleles,
        c_res,
        j_res,
        low_conf_js,
    ) = create_input_args(args, gene_types)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            stitched = st.stitch(
                input_args,
                tcr_dat,
                tcr_functionality,
                partial,
                codons,
                input_args['j_warning_threshold'],
                preferred_alleles,
                c_res,
                j_res,
                low_conf_js,
            )

    except ValueError:
        if cdr3[-1] == 'F':
            cdr3 = cdr3[:-1] + 'W'

        elif cdr3[-1] == 'W':
            cdr3 = cdr3[:-1] + 'F'

        args['cdr3'] = cdr3

        (
            input_args,
            tcr_dat,
            tcr_functionality,
            partial,
            codons,
            preferred_alleles,
            c_res,
            j_res,
            low_conf_js,
        ) = create_input_args(args, gene_types)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            stitched = st.stitch(
                input_args,
                tcr_dat,
                tcr_functionality,
                partial,
                codons,
                input_args['j_warning_threshold'],
                preferred_alleles,
                c_res,
                j_res,
                low_conf_js,
            )

    return fxn.translate_nt('N' * stitched['translation_offset'] + stitched['stitched_nt'])


def stcrpy_to_pandas(structure: Entity.Entity) -> pd.DataFrame:
    """Convert a stcrpy structure to a pandas dataframe."""
    records = [
        {
            'record_type': 'ATOM' if residue.id[0] == ' ' else 'HETATM',
            'atom_number': atom.serial_number,
            'atom_name': atom.id,
            'alt_loc': atom.altloc if atom.altloc != ' ' else None,
            'residue_name': residue.resname,
            'chain_id': chain.id,
            'residue_seq_id': residue.id[1],
            'residue_insert_code': residue.id[2] if residue.id[2] != ' ' else None,
            'pos_x': atom.coord[0],
            'pos_y': atom.coord[1],
            'pos_z': atom.coord[2],
            'occupancy': atom.occupancy,
            'b_factor': atom.bfactor,
            'element': atom.element,
            'charge': atom.get_charge(),
        }
        for chain in structure
        for residue in chain
        for atom in residue
    ]

    column_names = [
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

    return pd.DataFrame(records, columns=column_names)


def _merge_groups(groups: list[set[int]]) -> list[set[int]]:
    merged_groups = []

    for group in groups:
        merged = False

        for i, merged_group in enumerate(merged_groups):
            if len(group & merged_group) > 0:
                merged_groups[i] = group | merged_group
                merged = True
                break

        if not merged:
            merged_groups.append(group)

    return merged_groups


def find_common_groups(data: pd.DataFrame, columns: list[str]) -> np.ndarray:
    """Find common groups in a data frame based on input columns.

    Examples:
        >>> data
          col_a col_b  col_c
        0     a     l      1
        1     a     m      1
        2     b     l      2
        3     c     o      3
        4     c     p      0
        5     c     q      1
        6     d     q      5
        7     d     r      6
        8     e     s      1
        9     e     t     10
        >>> find_common_groups(data, ['col_a', 'col_b'])
        array([1, 1, 1, 2, 2, 2, 2, 2, 3, 3])

    Args:
        data: data frame containing the data to be grouped
        columns: columns to group based on shared values

    Returns:
        array with numbers corresponding to shared group ids

    Raises:
        ValueError: if no columns are input

    """
    if len(columns) == 0:
        msg = 'No columns input. Please input at least one column.'
        raise ValueError(msg)

    initial_groups = list(data.groupby(columns[0]))

    merge_matrix = np.zeros((len(initial_groups), len(initial_groups)), dtype=int)
    for i in range(len(initial_groups)):
        merge_matrix[i, i] = 1

    for i, (_, group_i) in enumerate(initial_groups):
        for j, (_, group_j) in enumerate(initial_groups[i + 1 :], i + 1):
            for column in columns[1:]:
                group_i_other = set(group_i[column].tolist())
                group_j_other = set(group_j[column].tolist())

                if len(group_i_other & group_j_other) > 0:
                    merge_matrix[i, j] = 1
                    break

    merge_matrix = np.maximum(merge_matrix, merge_matrix.T)

    group_indicies = np.arange(len(initial_groups))
    groups = [set(group_indicies[row > 0]) for row in merge_matrix]
    combined_groups = _merge_groups(groups)

    converted_groups = [
        {initial_groups[group_idx][0] for group_idx in grouped_groups} for grouped_groups in combined_groups
    ]
    group_map = {val: group_id for group_id, grouped_groups in enumerate(converted_groups, 1) for val in grouped_groups}

    return data[columns[0]].map(group_map).to_numpy()
