"""Select TCR-pMHC structures from STCRDab."""

import argparse
import logging
import os
import sys
import tempfile

import pandas as pd
import requests
from Bio.PDB import PDBIO, PDBParser, Structure
from sklearn.cluster import AgglomerativeClustering

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.comparisons import compute_structural_distances
from tcr_pmhc_predicted_structures.imgt_numbering import (
    IMGT_CDR,
    IMGT_CDR1,
    IMGT_CDR2,
    IMGT_CDR3,
    IMGT_FRAMEWORK_REGION,
    IMGT_MH1_ABD,
    IMGT_MH2_ABD,
    renumber_chain,
)
from tcr_pmhc_predicted_structures.missing_entities import (
    annotate_with_missing_entities,
    get_alignment,
    get_missing_atoms,
    get_missing_residues,
    screen_chain,
)
from tcr_pmhc_predicted_structures.missing_entities.fix import predict_missing_residues
from tcr_pmhc_predicted_structures.structure import (
    bio_to_pandas,
    crop_structure,
    extract_chains,
    get_header,
    get_sequence,
    remove_het_atoms,
    replace_chain,
)

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('stcrdab', help='path to the STCRDab')

parser.add_argument('--output', '-o', required=True, help='directory to output selected structures')
parser.add_argument('--output-summary-csv', required=True, help='path to output summary csv file of structures')

structure_type_group = parser.add_argument_group('Structure Types')
structure_type_group.add_argument(
    '--tcr-types',
    nargs='+',
    default=['abTCR'],
    help='TCR types allowed in dataset (abTCR and/or gdTCR) (Default: abTCR)',
)
structure_type_group.add_argument(
    '--mhc-types',
    nargs='+',
    default=['MH1'],
    help='MHC types allowed in dataset (CD1, GA, GB, MH1, MH2, MR1) (Default: MH1)',
)
structure_type_group.add_argument(
    '--antigen-types',
    nargs='+',
    default=['peptide'],
    help=('MHC types allowed in dataset (carbohydrate, Hapten, peptide, protein, etc) (Default: peptide)'),
)
structure_type_group.add_argument(
    '--crop-structures',
    action='store_true',
    help=('Crop TCR:pMHC structures to the TCR variable domain and MHC antigen binding domain.'),
)
structure_type_group.add_argument('--remove-het-atoms', action='store_true', help='Remove heteroatoms from structures')

quality_group = parser.add_argument_group('Quality Selection')
quality_group.add_argument(
    '--resolution-cutoff', type=float, default=3.50, help='maximum resolution allowed (Default: 3.50)'
)
quality_group.add_argument(
    '--remove-structures-missing-residues',
    action='store_true',
    help=(
        'Remove TCR:pMHC structures with missing residues in the TCR variable region or '
        'pMHC antigen binding domain (including peptide)'
    ),
)
quality_group.add_argument(
    '--fix-structures-missing-residues',
    action='store_true',
    help=(
        'Fix TCR:pMHC structures with missing residues in the TCR Fw region or '
        'MHC antigen binding domain rather than discarding (to be used with '
        '`--remove-structures-missing-residues)'
    ),
)
quality_group.add_argument(
    '--structural-similarity-cutoff',
    type=float,
    default=None,
    help='RMSD threshold for structures with the same CDR and peptide sequences (Default: None)',
)

add_logging_arguments(parser)

HISTO_FYI_BASE_URL = 'https://api.histo.fyi/v1/structures/{pdb}'


def fix_structure(
    structure: Structure.Structure,
    raw_structure: Structure.Structure,
    missing_residues: list[dict],
    missing_atoms: list[dict],
    chain_id: str,
) -> Structure.Structure:
    """Fix structure with missing residues or atoms on chain."""
    raw_structure_df = bio_to_pandas(raw_structure)
    raw_structure_df = annotate_with_missing_entities(raw_structure_df, missing_residues, missing_atoms)

    chain_df = raw_structure_df[raw_structure_df['chain_id'] == chain_id].query("record_type == 'ATOM'").copy()

    residue_numbering = chain_df[['residue_seq_id', 'residue_insert_code']].drop_duplicates()
    residue_numbering['res_num'] = range(1, len(residue_numbering) + 1)
    chain_df = chain_df.merge(residue_numbering, how='left', on=['residue_seq_id', 'residue_insert_code'])

    chain_missing_residues = chain_df.groupby('res_num').filter(lambda group: group['missing'].all())
    chain_missing_residues = chain_missing_residues.drop_duplicates('res_num')

    missing_residue_selection = chain_missing_residues.apply(
        lambda row: f"{row['res_num']}:{row['chain_id']}",
        axis=1,
    ).tolist()

    chain_missing_atoms = chain_df.groupby('res_num').filter(
        lambda group: group['missing'].any() and not group['missing'].all(),
    )
    chain_missing_atoms = chain_missing_atoms[chain_missing_atoms['missing']]
    missing_atom_selection = chain_missing_atoms.apply(
        lambda row: f"{row['atom_name']}:{row['res_num']}:{row['chain_id']}",
        axis=1,
    )

    alignment = get_alignment(structure, raw_structure, chain_id, missing_residues)

    new_chain = predict_missing_residues(
        alignment,
        structure,
        missing_residue_selection,
        missing_atom_selection,
        chain_id,
    )
    new_chain, _ = renumber_chain(new_chain)
    structure = replace_chain(structure, new_chain)

    return structure


def screen_for_missing_residues(
    df: pd.DataFrame, *, fix_structures: bool = False, fix_dir: str | None = None
) -> pd.DataFrame:
    """Remove entries missing residues in the TCR variable region or pMHC antigen binding domain."""

    def check_structure(  # noqa: PLR0911 TODO refactor this function
        pdb_id,
        alpha_chain_id,
        beta_chain_id,
        antigen_chain_id,
        mhc_chain1_id,
        mhc_chain2_id,
        mhc_type,
        imgt_file_path,
        raw_file_path,
    ):
        if mhc_type == 'MH1':
            mhc_chains = (mhc_chain1_id,)

        elif mhc_type == 'MH2':
            mhc_chains = (mhc_chain1_id, mhc_chain2_id)

        else:
            msg = f"Invalid MHC type {mhc_type}. MHC type must be 'MH1' or 'MH2'"
            raise ValueError(msg)

        logger.debug(
            'Checking %s, chains: %s, and MHC type %s',
            pdb_id,
            '-'.join([alpha_chain_id, beta_chain_id, antigen_chain_id, *mhc_chains]),
            mhc_type,
        )

        with open(raw_file_path, 'r') as fh:
            header = get_header(fh.read())

        pdb_parser = PDBParser(QUIET=True)
        raw_structure = pdb_parser.get_structure(pdb_id, raw_file_path)
        structure = pdb_parser.get_structure(pdb_id + '_imgt_numbered', imgt_file_path)

        missing_atoms = get_missing_atoms(header)
        missing_residues = get_missing_residues(header)

        fixed = False

        try:
            logger.debug('Checking for missing residues or atoms in antigen chain (chain %s)', antigen_chain_id)
            if antigen_chain_id in (
                {res['chain_id'] for res in missing_residues} | {res['chain_id'] for res in missing_atoms}
            ):
                logger.debug('Missing elements found')
                return 'missing'

            for tcr_chain_id in alpha_chain_id, beta_chain_id:
                logger.debug('Checking for missing residues or atoms in CDRs of chain %s', tcr_chain_id)
                if not screen_chain(
                    structure, raw_structure, tcr_chain_id, missing_residues, missing_atoms, IMGT_CDR, trim_ends=False
                ):
                    logger.debug('Missing elements found')
                    return 'missing'

                logger.debug('Checking for missing residues or atoms in framework region of chain %s', tcr_chain_id)
                if not screen_chain(
                    structure, raw_structure, tcr_chain_id, missing_residues, missing_atoms, IMGT_FRAMEWORK_REGION
                ):
                    logger.debug('Missing elements found')

                    if fix_structures:
                        logger.debug('Fixing chain')
                        structure = fix_structure(
                            structure,
                            raw_structure,
                            missing_residues,
                            missing_atoms,
                            tcr_chain_id,
                        )

                        fixed = True

                    else:
                        return 'missing'

            for mhc_chain_id in mhc_chains:
                logger.debug(
                    'Checking for missing residues or atoms in antigen binding domain of chain %s', mhc_chain_id
                )
                if not screen_chain(
                    structure,
                    raw_structure,
                    mhc_chain_id,
                    missing_residues,
                    missing_atoms,
                    IMGT_MH1_ABD if mhc_type == 'MH1' else IMGT_MH2_ABD,
                ):
                    logger.debug('Missing elements found')

                    if fix_structures:
                        logger.debug('Fixing chain')
                        structure = fix_structure(
                            structure,
                            raw_structure,
                            missing_residues,
                            missing_atoms,
                            mhc_chain_id,
                        )

                        fixed = True

                    else:
                        return 'missing'

            if fixed:
                io = PDBIO()
                io.set_structure(structure)
                io.save(os.path.join(fix_dir, os.path.basename(imgt_file_path)))

                return 'fixed'

        except KeyError:
            logger.warning(
                'Chain ID not found in raw structure of %s and chains %s. Structure will be removed.',
                pdb_id,
                '-'.join([alpha_chain_id, beta_chain_id, antigen_chain_id, *mhc_chains]),
            )
            return 'unknown'

        else:
            return 'complete'

    structure_status = df.apply(
        lambda row: check_structure(
            row.pdb,
            row.Achain,
            row.Bchain,
            row.antigen_chain,
            row.mhc_chain1,
            row.mhc_chain2,
            row.mhc_type,
            row.imgt_file_path,
            row.raw_file_path,
        ),
        axis=1,
    )

    return (
        df[(structure_status != 'missing') & (structure_status != 'unknown')].copy(),
        structure_status[(structure_status != 'missing') & (structure_status != 'unknown')].copy(),
    )


def add_cdr_sequences(pdb_id: str, alpha_chain_id: str, beta_chain_id: str, imgt_file_path: str) -> pd.Series:
    """Add CDR sequences from structures."""
    structure = PDBParser().get_structure(pdb_id, imgt_file_path)

    sequences = {
        f"cdr{cdr_num}_{chain_type.replace('_chain', '')}": get_sequence(structure, chain_id, numbering)
        for chain_type, chain_id in (('alpha_chain', alpha_chain_id), ('beta_chain', beta_chain_id))
        for cdr_num, numbering in ((1, IMGT_CDR1), (2, IMGT_CDR2), (3, IMGT_CDR3))
    }

    return pd.Series(sequences).sort_index()


def add_peptide_sequences(pdb_id: str, antigen_chain_id: str, imgt_file_path: str) -> str:
    """Add peptide sequences from structures."""
    structure = PDBParser().get_structure(pdb_id, imgt_file_path)

    return get_sequence(structure, antigen_chain_id)


def add_mhc_allele_information(pdb_id: str) -> tuple[str | None, str | None]:
    """Get MHC allele information from histo.fyi if available."""
    request = requests.get(HISTO_FYI_BASE_URL.format(pdb=pdb_id), timeout=20)
    content = request.json()

    if content['structure'] is None:
        return None, None

    alleles = tuple([content['structure']['allele'][chain].get('name', None) for chain in ('alpha', 'beta')])

    return alleles


def remove_similar_structures(df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Remove structures with the same CDR and peptide sequences within the RMSD threshold.

    The highest resolution structure will be kept.
    """
    output_dfs = []

    for (cdr_sequence, peptide_sequence, mhc_type), group in df.groupby(['collated_cdrs', 'peptide', 'mhc_type']):
        if len(group) == 1:
            output_dfs.append(group)
            continue

        logger.debug('Screening TCR: %s, peptide: %s, MHC: %s', cdr_sequence, peptide_sequence, mhc_type)

        sorted_group = group.sort_values(
            ['resolution', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2'],
        )

        pdb_parser = PDBParser(QUIET=True)

        structures = []
        chain_maps = []

        for _, row in sorted_group.iterrows():
            chain_map = {'alpha_chain': row.Achain, 'beta_chain': row.Bchain, 'antigen_chain': row.antigen_chain}

            if mhc_type == 'MH1':
                chain_map['mhc_chain1'] = row.mhc_chain1

            elif mhc_type == 'MH2':
                chain_map['mhc_chain1'] = row.mhc_chain1
                chain_map['mhc_chain2'] = row.mhc_chain2

            else:
                msg = f"Invalid MHC type {mhc_type}. MHC type must be 'MH1' or 'MH2'"
                raise ValueError(msg)

            logger.debug('Collecting PDB ID: %s and extracting chains %s', row.pdb, '-'.join(chain_map.values()))

            structure = pdb_parser.get_structure('', row.imgt_file_path)
            structure = extract_chains(structure, chain_map.values())

            structures.append(structure)
            chain_maps.append(chain_map)

        distance_matrix = compute_structural_distances(structures, chain_maps, mhc_type)
        clusters = (
            AgglomerativeClustering(
                metric='precomputed', distance_threshold=threshold, linkage='single', n_clusters=None
            )
            .fit(distance_matrix)
            .labels_
        )
        clusters = pd.Series(clusters, index=sorted_group.index)
        output_dfs.append(sorted_group[~clusters.duplicated()])

    return pd.concat(output_dfs)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    if args.fix_structures_missing_residues:
        fix_dir = tempfile.TemporaryDirectory()
        fix_dir_name = fix_dir.name

    else:
        fix_dir_name = None

    logger.info('Loading STCRDab summary')
    stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab, 'db_summary.dat'), delimiter='\t')
    logger.info('Number of structures: %d', len(stcrdab_summary))

    stcrdab_summary['imgt_file_path'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'imgt', pdb_id + '.pdb'),
    )
    stcrdab_summary['raw_file_path'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'raw', pdb_id + '.pdb'),
    )

    logger.info(
        'Selecting structures (TCRs: %s, MHCs: %s, antigens: %s)',
        ', '.join(args.tcr_types),
        ', '.join(args.mhc_types),
        ', '.join(args.antigen_types),
    )
    selected_structures = stcrdab_summary
    selected_structures = selected_structures.query('TCRtype in @args.tcr_types')
    selected_structures = selected_structures.query('mhc_type in @args.mhc_types')
    selected_structures = selected_structures.query('antigen_type in @args.antigen_types')
    selected_structures = selected_structures.copy()

    logger.info('Number of structures: %d', len(selected_structures))

    logger.info('Screening structures above %.2f Å resolution', args.resolution_cutoff)
    selected_structures['resolution'] = pd.to_numeric(selected_structures['resolution'], errors='coerce')
    selected_structures = selected_structures.query('resolution <= @args.resolution_cutoff')

    logger.info('Number of structures: %d', len(selected_structures))

    if args.remove_structures_missing_residues:
        logger.info('Removing structures missing residues')
        selected_structures, structure_status = screen_for_missing_residues(
            selected_structures,
            fix_structures=args.fix_structures_missing_residues,
            fix_dir=fix_dir_name,
        )

        logger.info('Number of structures: %d', len(selected_structures))

        if args.fix_structures_missing_residues:
            selected_structures.loc[structure_status == 'fixed', 'imgt_file_path'] = selected_structures[
                structure_status == 'fixed'
            ]['imgt_file_path'].map(
                lambda path: os.path.join(fix_dir_name, os.path.basename(path)),
            )

    logger.info('Getting sequence information...')
    cdr_sequences = selected_structures.apply(
        lambda row: add_cdr_sequences(row.pdb, row.Achain, row.Bchain, row.imgt_file_path),
        axis=1,
    )

    peptide_sequences = selected_structures.apply(
        lambda row: add_peptide_sequences(row.pdb, row.antigen_chain, row.imgt_file_path),
        axis=1,
    )
    peptide_sequences.name = 'peptide'

    selected_structures = pd.concat(
        [selected_structures, cdr_sequences, peptide_sequences.to_frame()],
        axis='columns',
    )

    selected_structures['collated_cdrs'] = (
        selected_structures['cdr1_alpha']
        + '-'
        + selected_structures['cdr2_alpha']
        + '-'
        + selected_structures['cdr3_alpha']
        + '-'
        + selected_structures['cdr1_beta']
        + '-'
        + selected_structures['cdr2_beta']
        + '-'
        + selected_structures['cdr3_beta']
    )

    logger.info('Getting available MHC information')
    selected_structures[['mhc1', 'mhc2']] = selected_structures['pdb'].map(add_mhc_allele_information).apply(pd.Series)
    selected_structures.loc[
        (selected_structures['mhc_type'] == 'MH1')
        & selected_structures['mhc2'].isna()
        & selected_structures['mhc1'].notna(),
        'mhc2',
    ] = 'B2M'

    if args.structural_similarity_cutoff:
        logger.info('Removing structures within %.2f Å RMSD', args.structural_similarity_cutoff)
        selected_structures = remove_similar_structures(selected_structures, args.structural_similarity_cutoff)

        logger.info('Number of structures: %d', len(selected_structures))

    selected_structures = selected_structures.fillna('')
    selected_structures = selected_structures.sort_values(
        [
            'pdb',
            'Achain',
            'Bchain',
            'antigen_chain',
            'mhc_chain1',
            'mhc_chain2',
        ]
    )
    selected_structures = selected_structures.reset_index()

    logger.info('Outputing structures...')
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    logger.debug('Formatting meta data')
    selected_structures['path'] = selected_structures.apply(
        lambda row: f'{row.pdb}_{row.Achain}{row.Bchain}{row.antigen_chain}{row.mhc_chain1}{row.mhc_chain2}.pdb',
        axis=1,
    )

    selected_structures = selected_structures.rename(
        {'alpha_subgroup': 'v_alpha', 'beta_subgroup': 'v_beta'},
        axis='columns',
    )
    selected_structures['species'] = (
        selected_structures.filter(regex='(?<!antigen_)organism')
        .fillna('')
        .apply(set, axis='columns')
        .map(list)
        .map(sorted)
        .map('/'.join)
        .str.strip('/')
        .str.replace('homo sapiens', 'Human')
        .str.replace('mus musculus', 'Mouse')
        .str.replace('triticum aestivum', 'Wheat')
        .str.replace('macaca mulatta', 'Rhesus Monkey')
        .str.replace('manduca sexta', 'Tobacco Hornworm')
        .str.replace('pan troglodytes', 'Chimpanzee')
        .str.replace('escherichia coli', 'E. coli')
    )

    selected_structures['mhc_chain2'] = selected_structures.apply(
        lambda row: '' if row.mhc_type == 'MH1' and args.crop_structures else row.mhc_chain2,
        axis=1,
    )

    output_columns = [
        'path',
        'pdb',
        'Achain',
        'Bchain',
        'antigen_chain',
        'mhc_chain1',
        'mhc_chain2',
        'mhc_type',
        'cdr1_alpha',
        'cdr2_alpha',
        'cdr3_alpha',
        'cdr1_beta',
        'cdr2_beta',
        'cdr3_beta',
        'v_alpha',
        'v_beta',
        'peptide',
        'mhc1',
        'mhc2',
        'species',
    ]

    logger.debug('Output summary file to %s', args.output_summary_csv)
    selected_structures[output_columns].to_csv(args.output_summary_csv, index=False)

    pdb_parser = PDBParser()
    for _, row in selected_structures.iterrows():
        structure = pdb_parser.get_structure(row.pdb, row.imgt_file_path)
        output_chains = [row.Achain, row.Bchain, row.antigen_chain, row.mhc_chain1, row.mhc_chain2]

        logger.debug('Outputting PDB structure %s...', row.path)

        if args.remove_het_atoms:
            logger.debug('Removing hetero atoms')
            structure = remove_het_atoms(structure)

        if args.crop_structures:
            logger.debug('Cropping structure')
            match row.mhc_type:
                case 'MH1':
                    structure = crop_structure(structure, (row.Achain, row.Bchain), (row.mhc_chain1,), row.mhc_type)

                    output_chains = [row.Achain, row.Bchain, row.antigen_chain, row.mhc_chain1]

                case 'MH2':
                    structure = crop_structure(
                        structure, (row.Achain, row.Bchain), (row.mhc_chain1, row.mhc_chain2), row.mhc_type
                    )

        structure = extract_chains(structure, output_chains)

        io = PDBIO()
        io.set_structure(structure)
        io.save(os.path.join(args.output, row.path))

    if args.fix_structures_missing_residues:
        fix_dir.cleanup()


if __name__ == '__main__':
    main()
