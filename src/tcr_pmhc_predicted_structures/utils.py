"""Utility functions for processing TCR:pMHC data."""

import logging
import re
import warnings
from functools import lru_cache

import pandas as pd
from Bio.PDB import Entity
from Stitchr import stitchr as st
from Stitchr import stitchrfunctions as fxn

logger = logging.getLogger(__name__)

IMGT_VARIABLE_DOMAIN: set[int] = set(range(1, 128 + 1))
'''Variable domain range for IMGT numbered TCR structures.'''

IMGT_CDR1: set[int] = set(range(27, 38 + 1))
'''IMGT residue numbers corresponding to CDR 1 domains.'''
IMGT_CDR2: set[int] = set(range(56, 65 + 1))
'''IMGT residue numbers corresponding to CDR 2 domains.'''
IMGT_CDR3: set[int] = set(range(105, 117 + 1))
'''IMGT residue numbers corresponding to CDR 3 domains.'''
IMGT_CDR: set[int] = IMGT_CDR1.union(IMGT_CDR2).union(IMGT_CDR3)
'''IMGT residue numbers corresponding to all CDR domains.'''

IMGT_MH1_ABD: set[int] = set(range(1, 92)) | set(range(1001, 1092))
'''IMGT ranges of the antigen binding domain of MHC class I molecules.'''

IMGT_MH2_ABD: set[int] = set(range(1, 92))
'''IMGT ranges of the antigen binding domain of MHC class II molecules.'''

MHC_I_IMGT_BETA_HELIX_START: int = 1000
'''Start of the beta helix of class I molecules.'''

IMGT_MHC_TCR_CONTACT_POSITIONS = {
    'alpha': ['62', '63', '65', '66', '68', '69', '70', '72', '73', '75', '76', '79'],
    'beta': ['58', '61A', '62', '63', '65', '66', '67', '69', '70', '72', '72A', '73', '76', '77'],
}
'''IMGT Residues contacted by TCRs on the surface of MHC molecules.'''


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
