"""Convert sequence data into different formats for structure prediction methods."""

import argparse
import csv
import glob
import json
import logging
import os
import re
import sys

import pandas as pd

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_CDR1, IMGT_CDR2, IMGT_CDR3, MHC_I_IMGT_BETA_HELIX_START

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
parser.add_argument('--chai-restraints', help='restraints file for chai restraints')

paths_group = parser.add_argument_group('Paths')
paths_group.add_argument('--sequence-key-path', help='path to sequence key file')
paths_group.add_argument(
    '--relative-paths',
    action='store_true',
    help='use relative path names for MSA and template paths instead of full path names',
)
paths_group.add_argument('--unpaired-msa-paths', nargs='+', help='paths to the unpaired MSA files')
paths_group.add_argument('--paired-msa-paths', nargs='+', help='paths to the paired MSA files')
paths_group.add_argument('--template-paths', nargs='+', help='path to the templates')

add_logging_arguments(parser)

CHAI_CHAIN_ORDER = ['mhc_chain_1', 'mhc_chain_2', 'peptide', 'tcr_alpha_chain', 'tcr_beta_chain']


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


def format_sequences_for_af3(
    sequences: list[dict],
    output_dir: str,
    model_seeds: list[int],
    sequence_key: dict[str, str] | None = None,
    unpaired_msa_paths: list[str] | None = None,
    paired_msa_paths: list[str] | None = None,
    templates: dict[str, list[dict[str, str]]] | None = None,
) -> None:
    if unpaired_msa_paths is None:
        unpaired_msa_paths = []

    if paired_msa_paths is None:
        paired_msa_paths = []

    if templates is None:
        templates = []

    for sequence_info in sequences:
        logger.debug('Writing json for %s', sequence_info['id'])
        output_name = os.path.join(output_dir, f"{sequence_info['id']}_{'-'.join(model_seeds)}.json")

        output = {
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
            'version': 2,
        }

        if sequence_key:
            new_sequence_entries = []
            for sequence_entry in output['sequences']:
                new_sequence_entry = {**sequence_entry}

                if sequence_entry['protein']['id'] == 'C':
                    msa = f">query\n{sequence_entry['protein']['sequence']}\n"
                    new_sequence_entry['protein']['unpairedMsa'] = msa
                    new_sequence_entry['protein']['pairedMsa'] = msa
                    new_sequence_entry['protein']['templates'] = []

                else:
                    try:
                        key = sequence_key[sequence_entry['protein']['sequence']]

                    except KeyError:
                        logger.error('No entry found for sequence: %s', sequence_entry['protein']['sequence'])  # noqa: TRY400

                    else:
                        relevant_unpaired_msas = [msa_path for msa_path in unpaired_msa_paths if key in msa_path]

                        if len(relevant_unpaired_msas) < 1:
                            logger.error(
                                'Missing unpaired MSAs for sequence: %s (key: %s)',
                                sequence_entry['protein']['sequence'],
                                key,
                            )

                        else:
                            new_sequence_entry['protein']['unpairedMsaPath'] = relevant_unpaired_msas[0]

                        relevant_paired_msas = [msa_path for msa_path in paired_msa_paths if key in msa_path]

                        if len(relevant_unpaired_msas) < 1:
                            logger.warning(
                                'Missing paired MSAs for sequence: %s (key: %s)',
                                sequence_entry['protein']['sequence'],
                                key,
                            )

                        else:
                            new_sequence_entry['protein']['pairedMsaPath'] = relevant_paired_msas[0]

                        try:
                            new_sequence_entry['protein']['templates'] = templates[key]

                        except KeyError:
                            logger.warning(
                                'Missing paired MSAs for sequence: %s (key: %s)',
                                sequence_entry['protein']['sequence'],
                                key,
                            )

                new_sequence_entries.append(new_sequence_entry)

            output['sequences'] = new_sequence_entries

        with open(output_name, 'w') as output_fh:
            logger.debug('Writing to %s', output_name)
            json.dump(output, output_fh)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading sequences from %s', args.input)
    with open(args.input, 'r') as fh:
        sequences = list(csv.DictReader(fh))

    if args.sequence_key_path:
        logger.info('Loading sequence keys from %s', args.sequence_key_path)
        with open(args.sequence_key_path, 'r') as fh:
            sequence_key = {sequence: key for key, sequence in (line.split() for line in fh if line)}

        unpaired_msa_paths = args.unpaired_msa_paths
        if args.relative_paths and unpaired_msa_paths:
            dirname = os.path.dirname(unpaired_msa_paths[0])
            unpaired_msa_paths = [path.replace(dirname, '').strip('/') for path in unpaired_msa_paths]

        paired_msa_paths = args.paired_msa_paths
        if args.relative_paths and paired_msa_paths:
            dirname = os.path.dirname(paired_msa_paths[0])
            paired_msa_paths = [path.replace(dirname, '').strip('/') for path in paired_msa_paths]

        if args.template_paths:
            templates = {}
            for path in args.template_paths:
                name = os.path.basename(path)

                indices_paths = sorted(glob.glob(os.path.join(path, '*_indices.json')))
                structure_paths = sorted(glob.glob(os.path.join(path, '*.cif')))

                entry_templates = []
                for indices_path, structure_path in zip(indices_paths, structure_paths, strict=True):
                    with open(indices_path, 'r') as fh:
                        entry_template = json.load(fh)

                    if args.relative_paths:
                        dirname = os.path.dirname(path)
                        entry_template['mmcifPath'] = structure_path.replace(dirname, '').strip('/')

                    else:
                        entry_template['mmcifPath'] = structure_path

                    entry_templates.append(entry_template)

                templates[name] = entry_templates

        else:
            templates = None

    else:
        sequence_key = None
        unpaired_msa_paths = None
        paired_msa_paths = None
        templates = None

    if not os.path.exists(args.output):
        logger.debug('Creating output dir: %s', args.output)
        os.mkdir(args.output)

    match args.format:
        case 'alphafold3':
            logger.info('Formatting sequences for Alphafold3 and outputting to %s', args.output)
            format_sequences_for_af3(
                sequences,
                args.output,
                args.alphafold3_model_seeds,
                sequence_key,
                unpaired_msa_paths,
                paired_msa_paths,
                templates,
            )

        case 'chai':
            logger.info('Formatting sequences for chai and outputting to %s', args.output)

            if args.chai_restraints:
                logger.info('Loading restraints from %s', args.chai_restraints)
                restraints = pd.read_csv(args.chai_restraints)

            else:
                restraints = None

            format_sequences_for_chai(sequences, args.output, restraints)


if __name__ == '__main__':
    main()
