"""Constants and functions for annotating sequences using IMGT numbering."""

import logging

from Bio.Align import PairwiseAligner
from Bio.PDB import Chain
from Bio.SeqUtils import IUPACData

logger = logging.getLogger(__name__)

try:
    import anarci

except ImportError:
    logger.warning('Anarci is unavailable and is required by some functions')

IMGT_CDR1: set[int] = set(range(27, 38 + 1))
'''IMGT residue numbers corresponding to CDR 1 domains.'''
IMGT_CDR2: set[int] = set(range(56, 65 + 1))
'''IMGT residue numbers corresponding to CDR 2 domains.'''
IMGT_CDR3: set[int] = set(range(105, 117 + 1))
'''IMGT residue numbers corresponding to CDR 3 domains.'''
IMGT_CDR: set[int] = IMGT_CDR1.union(IMGT_CDR2).union(IMGT_CDR3)
'''IMGT residue numbers corresponding to all CDR domains.'''

IMGT_VARIABLE_DOMAIN: set[int] = set(range(1, 128 + 1))
'''Variable domain range for IMGT numbered TCR structures.'''

IMGT_FRAMEWORK_REGION: set[int] = IMGT_VARIABLE_DOMAIN - IMGT_CDR
'''Framework (Fw) region range for IMGT numbered TCR structures.'''

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

AnarciNumbering = list[tuple[tuple[int, str], str]]

NUMBERING_OFFSET = 1000
'''Offset needed to number things like scTCRs where two numberings are present and they may clash.'''


def assign_cdr_number(seq_id: int) -> int | None:
    """Assign CDR number for a sequence ID or return None if the sequence ID if the ID is not in CDR range."""
    if seq_id in IMGT_CDR1:
        return 1

    if seq_id in IMGT_CDR2:
        return 2

    if seq_id in IMGT_CDR3:
        return 3

    return None


def _get_numberings(sequence: str) -> AnarciNumbering:
    seq_position = 0
    seq_length = len(sequence)
    finished_numbering = False

    numberings = []
    while not finished_numbering:
        numbering, chain_type = anarci.number(sequence[seq_position:])

        if not numbering:
            finished_numbering = True
            break

        if len(numbering) == 2:  # noqa: PLR2004 allow magic number here
            numbering = numbering[0]
            chain_type = chain_type[0]

            logger.warning(
                ('Multiple possible chain annotations found for sequence %s. Defaulting to first (%sCHAIN)'),
                sequence,
                chain_type,
            )

        clean_numbering = [(_, aa) for _, aa in numbering if aa != '-']

        # Needed for scTCRs
        clean_seq_ids = {seq_id for seq_id, _ in clean_numbering}
        numberings_seq_ids = {seq_id for seq_id, _ in numberings}

        if len(clean_seq_ids & numberings_seq_ids) > 0:
            clean_numbering = [
                ((res_seq_id + NUMBERING_OFFSET, insert_code), res_olc)
                for (res_seq_id, insert_code), res_olc in clean_numbering
            ]

        numberings.extend(clean_numbering)

        clean_sequence = ''.join([olc for _, olc in clean_numbering])
        seq_position += sequence.index(clean_sequence) + len(clean_sequence)

        if seq_position >= seq_length:
            finished_numbering = True

    return numberings


class AlignmentError(Exception):
    """Raised when there is an error aligning two sequences."""

    def __init__(self, sequence1: str, sequence2: str, details: str | None = None) -> None:
        message = f'Could not align sequences: {sequence1} and {sequence2}'

        if details:
            message += '. ' + details

        super().__init__(message)


def _align_numbering(sequence: str, numbering: AnarciNumbering) -> AnarciNumbering:
    numbering_sequence = ''.join([olc for _, olc in numbering])

    aligner = PairwiseAligner(match_score=1, open_gap_score=-1, end_gap_score=0.001)
    alignments = aligner.align(sequence, numbering_sequence)
    alignment = alignments[0]

    numbered_index = 0

    seq_id = (1, ' ')
    start_flag = True
    numbered_positions = {seq_id for seq_id, _ in numbering}

    full_numbering = []
    for input_index, (numbered_residue, input_residue) in enumerate(zip(alignment[1], alignment[0], strict=True)):
        if numbered_residue == input_residue:
            seq_id = numbering[numbered_index][0]

            numbered_index += 1
            start_flag = False

        elif numbered_residue == '-':
            propposed_seq_id = seq_id if start_flag else (seq_id[0] + 1, ' ')

            if propposed_seq_id in numbered_positions:
                for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    propposed_seq_id = (seq_id[0], letter)

                    if propposed_seq_id not in numbered_positions:
                        break

                else:
                    raise AlignmentError(
                        alignment[1],
                        alignment[0],
                        f'Not enough insertion codes to cover gap in sequences at position {seq_id[0]}.',
                    )

            seq_id = propposed_seq_id
            start_flag = False

        else:
            raise AlignmentError(alignment[1], alignment[0])

        full_numbering.append((seq_id, sequence[input_index]))
        numbered_positions.add(seq_id)

    return full_numbering


def number(sequence: str) -> AnarciNumbering:
    """Numbers sequence with IMGT convention using ANARCI."""
    numbering = _get_numberings(sequence)

    if not numbering:
        logger.warning('Sequence not identified by ANARCI: %s. Numbering sequence sequentially', {sequence})
        return [((seq_id, ' '), olc) for seq_id, olc in enumerate(sequence, 1)]

    return _align_numbering(sequence, numbering)


def renumber_chain(chain: Chain.Chain) -> tuple[Chain.Chain, str]:
    """Renumber a chain following IMGT conventions.

    Args:
        chain: input chain

    Returns:
        tuple with the new chain with imgt numbering and the type of chain identified by anarci

    Raises:
        ValueError: if chain cannot be identified by ANARCI

    """
    residues = list(chain.get_residues())

    sequence = ''.join(
        [
            IUPACData.protein_letters_3to1[res.get_resname().title()]
            for res in residues
            if res.id[0] == ' ' and res.get_resname().title() in IUPACData.protein_letters_3to1
        ]
    )

    numbering, chain_type = anarci.number(sequence)

    if not numbering:
        msg = f'Chain ID {chain.id} not identified by ANARCI'
        raise ValueError(msg)

    if len(numbering) == 2:  # noqa: PLR2004
        numbering = numbering[0]
        chain_type = chain_type[0]

        logger.warning(
            ('Multiple possible chain annotations found for chain id %s. Defaulting to first: %sCHAIN'),
            chain.id,
            chain_type,
        )

    numbering = [(seq_id, insert_code) for (seq_id, insert_code), res_name in numbering if res_name != '-']

    num_residues_not_numbered = len(residues) - len(numbering)

    next_seq_id = numbering[-1][0] + 1

    for _ in range(num_residues_not_numbered):
        numbering.append((next_seq_id, ' '))
        next_seq_id += 1

    output_chain = Chain.Chain(chain.id)
    for (seq_id, insert_code), res in zip(numbering, residues, strict=False):
        new_res = res.copy()
        new_res.id = (new_res.id[0], seq_id, insert_code)
        output_chain.add(new_res)

    return output_chain, chain_type
