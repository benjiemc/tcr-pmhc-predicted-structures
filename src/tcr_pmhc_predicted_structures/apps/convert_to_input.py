"""Convert sequence data into different formats for structure prediction methods."""

import argparse
import csv
import json
import logging
import os
import re
import sys

import pandas as pd
import tidytcells

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3, MHC_I_IMGT_BETA_HELIX_START
from tcr_pmhc_predicted_structures.utils import (
    assign_species,
    mhc_slug_to_code,
    stitch_sequence,
)

logger = logging.getLogger()

try:
    import anarci

except ImportError:
    logger.warning('Anarci is unavailable and is required by some functions')

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('input', help='path to input csv')
parser.add_argument('--output', '-o', required=True, help='')
parser.add_argument('--format', required=True, choices=['alphafold3', 'chai'], help='Model to generate input files for')
parser.add_argument(
    '--alphafold3-model-seeds', nargs='+', type=int, default=[1], help='Seeds for alphafold3 (Default: [1])'
)
parser.add_argument(
    '--mhc-sequences',
    nargs='+',
    required=True,
    help='paths to mhc sequence information to create pseudo sequences.',
)
parser.add_argument('--chai-restraints', help='restraints file for chai restraints')

add_logging_arguments(parser)

B2M_SEQUENCES = {
    'Human': 'MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPADIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM',
    'Mouse': (
        'MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEP'
        'KTVYWDRDM'
    ),
}
'''B2M sequences downloaded from UniProt.'''

CHAI_CHAIN_ORDER = ['mhc_chain_1', 'mhc_chain_2', 'peptide', 'tcr_alpha_chain', 'tcr_beta_chain']


def create_full_length_sequences(sequences: list[dict], mhc_sequences: pd.DataFrame) -> list[dict]:
    """Create full length TCR:pMHC sequences from V and J genes, CDR3s, MHC allele, and peptide sequence."""
    sequences_full = []

    for sequence in sequences:
        logger.debug('Working on %s', sequence['index'])
        sequence['v_alpha'] = sequence['v_alpha'].split(';')[0]
        sequence['j_alpha'] = sequence['j_alpha'].split(';')[0]
        sequence['v_beta'] = sequence['v_beta'].split(';')[0]
        sequence['j_beta'] = sequence['j_beta'].split(';')[0]

        logger.debug('Stitching TCR alpha chain')
        for species in sequence['species'].split('/'):
            tcr_alpha_seq = stitch_sequence(sequence['v_alpha'], sequence['j_alpha'], sequence['cdr3_alpha'], species)

            if tcr_alpha_seq is not None:
                break

        else:
            logger.warning(
                'Unable to stitch: %s-%s-%s, %s. Skipping...',
                sequence['v_alpha'],
                sequence['cdr3_alpha'],
                sequence['j_alpha'],
                sequence['tcr_species'],
            )
            continue

        logger.debug('Stitching TCR beta chain')
        for species in sequence['species'].split('/'):
            tcr_beta_seq = stitch_sequence(sequence['v_beta'], sequence['j_beta'], sequence['cdr3_beta'], species)

            if tcr_beta_seq is not None:
                break

        else:
            logger.warning(
                'Unable to stitch: %s-%s-%s, %s. Skipping...',
                sequence['v_beta'],
                sequence['cdr3_beta'],
                sequence['j_beta'],
                sequence['tcr_species'],
            )
            continue

        logger.debug('Looking up MHC chain 1 sequence')
        mhc_1_seq = mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc1']].iloc[0]['canonical_sequence']
        mhc_species = mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc1']].iloc[0]['species']

        logger.debug('Looking up MHC chain 2 sequence')
        mhc_2_seq = (
            B2M_SEQUENCES[mhc_species]
            if sequence['mhc2'] == 'B2M'
            else mhc_sequences[mhc_sequences['mhc_code'] == sequence['mhc2']].iloc[0]['canonical_sequence']
        )

        sequences_full.append(
            {
                'id': sequence['index'],
                'tcr_alpha_chain': tcr_alpha_seq,
                'tcr_beta_chain': tcr_beta_seq,
                'peptide': sequence['peptide'],
                'mhc_chain_1': mhc_1_seq,
                'mhc_chain_2': mhc_2_seq,
                'mhc_type': sequence['mhc_type'],
            }
        )

    return sequences_full


def create_restraints(sequence_info: dict, restraints: pd.DataFrame) -> pd.DataFrame:
    """Create restraints for a specific sequence.

    Args:
        sequence_info: information for the sequence
        restraints: general restraints for TCR:pMHCs

    Returns:
        restraints that are specific to the input sequence

    """
    chai_restraints = pd.DataFrame(
        {
            'restraint_id': restraints['cdr_anchor_tcr'] + '_' + restraints['mhc_tcr_contact_mhc'],
            'connection_type': ['contact'] * len(restraints),
            'confidence': [1.0] * len(restraints),
            'min_distance_angstrom': [0.0] * len(restraints),
            'max_distance_angstrom': restraints['distance'],
            'comment': ['blank'] * len(restraints),
        }
    )

    chai_restraints['chainA'] = [
        'D' if 'alpha' in cdr_location else 'E' for cdr_location in restraints['cdr_anchor_tcr']
    ]
    chai_restraints['chainB'] = [
        'A' if 'alpha' in mhc_tcr_contact_mhc else 'B' for mhc_tcr_contact_mhc in restraints['mhc_tcr_contact_mhc']
    ]

    if sequence_info['mhc_type'] == 'MH1':
        numbered_sequences = {
            name: anarci.number(sequence_info[name])[0]
            for name in CHAI_CHAIN_ORDER
            if name not in {'peptide', 'mhc_chain_2'}
        }

    else:
        numbered_sequences = {
            name: anarci.number(sequence_info[name])[0] for name in CHAI_CHAIN_ORDER if name not in {'peptide'}
        }

    numbered_sequences = {
        name: [
            ((seq_id, insert_code), res_olc) for (seq_id, insert_code), res_olc in numbered_sequence if res_olc != '-'
        ]
        for name, numbered_sequence in numbered_sequences.items()
    }

    res_idx_a = []
    for cdr_location in restraints['cdr_anchor_tcr']:
        cdr_num, chain_type, anchor = re.search('CDR([1-3])(alpha|beta)_(start|end)', cdr_location).group(1, 2, 3)
        cdr_imgt_range = (IMGT_CDR1, IMGT_CDR2, IMGT_CDR3)[int(cdr_num) - 1]

        match anchor:
            case 'start':
                prev_residue = {}
                for res_pos, ((seq_id, _), res_olc) in enumerate(numbered_sequences[f'tcr_{chain_type}_chain'], 1):
                    if seq_id in cdr_imgt_range:
                        res_idx_a.append(f"{prev_residue['res_olc']}{prev_residue['res_pos']}")
                        break

                    prev_residue['res_olc'] = res_olc
                    prev_residue['res_pos'] = res_pos

            case 'end':
                before_cdr = True
                in_cdr = False

                for res_pos, ((seq_id, _), res_olc) in enumerate(numbered_sequences[f'tcr_{chain_type}_chain'], 1):
                    if seq_id in cdr_imgt_range and before_cdr:
                        before_cdr = False
                        in_cdr = True

                    elif seq_id not in cdr_imgt_range and in_cdr:
                        res_idx_a.append(f'{res_olc}{res_pos}')
                        break

    chai_restraints['res_idxA'] = res_idx_a

    res_idx_b = []
    for mhc_tcr_contact in restraints['mhc_tcr_contact_mhc']:
        helix, imgt_seq_id, imgt_insert_code = re.search(r'(alpha|beta)_(\d+)([A-Za-z])?', mhc_tcr_contact).group(
            1, 2, 3
        )

        imgt_seq_id = int(imgt_seq_id)
        imgt_insert_code = imgt_insert_code if imgt_insert_code else ' '

        if sequence_info['mhc_type'] == 'MH1' and helix == 'beta':
            chain_name = 'mhc_chain_1'
            imgt_seq_id += MHC_I_IMGT_BETA_HELIX_START

        elif sequence_info['mhc_type'] == 'MH2' and helix == 'beta':
            chain_name = 'mhc_chain_2'

        else:
            chain_name = 'mhc_chain_1'

        for res_pos, ((seq_id, insert_code), res_olc) in enumerate(numbered_sequences[chain_name], 1):
            if imgt_seq_id == seq_id and imgt_insert_code == insert_code:
                res_idx_b.append(f'{res_olc}{res_pos}')
                break

    chai_restraints['res_idxB'] = res_idx_b

    return chai_restraints[
        [
            'restraint_id',
            'chainA',
            'res_idxA',
            'chainB',
            'res_idxB',
            'connection_type',
            'confidence',
            'min_distance_angstrom',
            'max_distance_angstrom',
            'comment',
        ]
    ]


def format_sequences_for_chai(sequences: list[dict], output_dir: str, restraints: pd.DataFrame | None = None) -> None:
    """Create the input files for chai.

    Args:
        sequences: list of dictionaries containing sequences for the tcr and pmhc chains
        output_dir: path to output fasta and constraint files to
        restraints: general restraints for TCR:pMHCs

    """
    for sequence_info in sequences:
        if restraints is not None:
            logger.debug('Creating restraints for %s', sequence_info['id'])
            sequence_restraints = create_restraints(sequence_info, restraints)

            restraints_path = os.path.join(output_dir, sequence_info['id'] + '_restraints.csv')
            logger.debug('Writing restraints CSV to %s', restraints_path)
            sequence_restraints.to_csv(restraints_path, index=False)

        logger.debug('Writing fasta for %s', sequence_info['id'])
        with open(os.path.join(output_dir, sequence_info['id'] + '.fasta'), 'w') as fh:
            for name in CHAI_CHAIN_ORDER:
                fh.write(f'>protein|name={name}\n')
                fh.write(sequence_info[name])
                fh.write('\n')


def format_sequences_for_af3(sequences: list[dict], output_dir: str, model_seeds: list[int]) -> None:
    for sequence_info in sequences:
        logger.debug('Writing json for %s', sequence_info['id'])
        output_name = os.path.join(output_dir, f"{sequence_info['id']}.json")
        with open(output_name, 'w') as output_fh:
            logger.debug('Writing to %s', output_name)
            json.dump(
                {
                    'name': sequence_info['id'],
                    'modelSeeds': model_seeds,
                    'sequences': [
                        {'protein': {'id': 'D', 'sequence': sequence_info['tcr_alpha_chain']}},
                        {'protein': {'id': 'E', 'sequence': sequence_info['tcr_beta_chain']}},
                        {'protein': {'id': 'C', 'sequence': sequence_info['peptide']}},
                        {'protein': {'id': 'A', 'sequence': sequence_info['mhc_chain_1']}},
                        {'protein': {'id': 'B', 'sequence': sequence_info['mhc_chain_2']}},
                    ],
                    'dialect': 'alphafold3',
                    'version': 1,
                },
                output_fh,
            )


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading MHC canonical sequences')
    mhc_sequences = pd.concat([pd.read_json(path, orient='index') for path in args.mhc_sequences])

    mhc_sequences['mhc_code'] = mhc_sequences.index.map(mhc_slug_to_code)
    mhc_sequences['species'] = mhc_sequences['mhc_code'].map(assign_species)

    mhc_sequences.loc[mhc_sequences['species'] == 'Human', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Human', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='homosapiens')

    mhc_sequences.loc[mhc_sequences['species'] == 'Mouse', 'mhc_code'] = mhc_sequences.loc[
        mhc_sequences['species'] == 'Mouse', 'mhc_code'
    ].apply(tidytcells.mh.standardize, species='musmusculus')

    logger.info('Loading sequences from %s', args.input)
    with open(args.input, 'r') as fh:
        sequences = list(csv.DictReader(fh))

    logger.info('Extending sequences')
    full_sequences = create_full_length_sequences(sequences, mhc_sequences)

    if not os.path.exists(args.output):
        logger.debug('Creating output dir: %s', args.output)
        os.mkdir(args.output)

    match args.format:
        case 'alphafold3':
            logger.info('Formatting sequences for Alphafold3 and outputting to %s', args.output)
            format_sequences_for_af3(full_sequences, args.output, args.alphafold3_model_seeds)

        case 'chai':
            logger.info('Formatting sequences for chai and outputting to %s', args.output)

            if args.chai_restraints:
                logger.info('Loading restraints from %s', args.chai_restraints)
                restraints = pd.read_csv(args.chai_restraints)

            else:
                restraints = None

            format_sequences_for_chai(full_sequences, args.output, restraints)


if __name__ == '__main__':
    main()
