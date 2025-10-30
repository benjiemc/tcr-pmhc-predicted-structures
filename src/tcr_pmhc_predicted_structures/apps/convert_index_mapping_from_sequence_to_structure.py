"""Convert the index mapping based on PDB Seqres records to atomic records."""

import argparse
import json
import logging
import sys

from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser
from Bio.SeqIO.PdbIO import CifSeqresIterator

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('structure', help='path to mmcif structure')
parser.add_argument('indices', help='path to indices to convert in AF3 json format (conversion is done for template)')
parser.add_argument('--output', '-o', required=True, help='path to output converted indices to (json format)')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading indices from %s', args.indices)
    with open(args.indices, 'r') as fh:
        indices = json.load(fh)

    logger.info('Loading sequence from Atom coordinates of %s', args.structure)
    structure = MMCIFParser(QUIET=True).get_structure('', args.structure)

    models = list(structure.get_models())
    if len(models) > 1:
        msg = 'Multiple models found, please ensure to provide a structure with only one model'
        raise ValueError(msg)

    chains = list(structure[0].get_chains())
    if len(chains) > 1:
        msg = 'Multiple chains found, please ensure to provide a structure with only one chain'
        raise ValueError(msg)

    structure_sequence = ''.join([protein_letters_3to1_extended[res.get_resname()] for res in chains[0]])
    logger.debug('Structure sequence: %s', structure_sequence)

    logger.info('Loading sequence from SeqRes of %s', args.structure)
    with open(args.structure, 'r') as fh:
        records = list(CifSeqresIterator(fh))
        if len(records) > 1:
            msg = 'Multiple SeqRes chains found, please ensure to provide a structure with only one SeqRes entry'
            raise ValueError(msg)

        seqres_sequence = str(records[0].seq)

    logger.debug('SeqRes sequence: %s', seqres_sequence)

    logger.info('Aligning sequences')

    aligner = PairwiseAligner(match_score=1, open_gap_score=-1, end_gap_score=0.001)
    alignments = aligner.align(structure_sequence, seqres_sequence)
    alignment = alignments[0]

    logger.debug('Alignement:\n\n%s', str(alignment))

    logger.info('Creating map between PDBSeqres indices and atomic indices')
    seqres_to_structure = {
        int(idx_seqres): int(idx_structure)
        for idx_seqres, idx_structure in zip(alignment.indices[1], alignment.indices[0], strict=True)
        if idx_seqres >= 0 and idx_structure >= 0
    }

    logger.debug(
        'Map between SeqRes and structure:\n\nSeqRes-Structure\n---------------\n%s',
        '\n'.join('-'.join(str(idx) for idx in idx_pair) for idx_pair in sorted(seqres_to_structure.items())),
    )

    logger.info('Converting SeqRes indices to structure indices')
    converted_indices = {
        query_idx: seqres_to_structure[template_idx]
        for query_idx, template_idx in zip(indices['queryIndices'], indices['templateIndices'], strict=True)
        if template_idx in seqres_to_structure
    }

    logger.info('Outputting to %s', args.output)
    with open(args.output, 'w') as fh:
        json.dump(
            {
                'queryIndices': list(converted_indices.keys()),
                'templateIndices': list(converted_indices.values()),
            },
            fh,
        )


if __name__ == '__main__':
    main()
