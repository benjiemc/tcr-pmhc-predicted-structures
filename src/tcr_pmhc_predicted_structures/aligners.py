"""Functions for aligning data."""

import numpy as np


def align_sequences(
    seq1: str, seq2: str, match_score: float = 1.0, mismatch_score: float = -1.0, indel_score: float = -1.0
) -> tuple:
    """Align two sequences using the Needleman-Wunch algorithm.

    See here https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm for more details.

    Args:
        seq1: first string to compare
        seq2: second string to compare
        match_score: score awarded for a match (Default: 1.0)
        mismatch_score: score awarded for a mismtach (Default: -1.0)
        indel_score: score awarded for an insertion or a deletion (Default: -1.0)

    Returns:
        tuple containing the alignment (represented as a list with tuple for each pair) followed by the score given to
        the alignment.

    """
    num_cols = len(seq1) + 1
    num_rows = len(seq2) + 1

    # Initialize the matrix
    alignment_matrix = np.empty((num_rows, num_cols))
    alignment_matrix[:] = np.nan

    alignment_matrix[0, 0] = 0

    alignment_matrix[1:, 0] = np.cumsum(np.repeat(indel_score, num_rows - 1))
    alignment_matrix[0, 1:] = np.cumsum(np.repeat(indel_score, num_cols - 1))

    direction_matrix = np.empty((num_rows, num_cols), dtype=str)

    direction_matrix[1:, 0] = 't'
    direction_matrix[0, 1:] = 'l'

    # Fill the matrix
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            top_score = alignment_matrix[i - 1, j] + indel_score
            left_score = alignment_matrix[i, j - 1] + indel_score
            diagonal_score = alignment_matrix[i - 1, j - 1] + (
                match_score if seq1[j - 1] == seq2[i - 1] else mismatch_score
            )
            scores = np.array((top_score, left_score, diagonal_score))
            score_taken = np.argmax(scores)

            alignment_matrix[i, j] = np.max(scores)
            direction_matrix[i, j] = ['t', 'l', 'd'][score_taken]

    # Trace optimal path
    reverse_alignment = []

    row_pos = num_rows - 1
    col_pos = num_cols - 1

    while (row_pos, col_pos) != (0, 0):
        if direction_matrix[row_pos, col_pos] == 't':
            reverse_alignment.append(('-', seq2[row_pos - 1]))

            row_pos -= 1

        elif direction_matrix[row_pos, col_pos] == 'l':
            reverse_alignment.append((seq1[col_pos - 1], '-'))

            col_pos -= 1

        else:
            reverse_alignment.append((seq1[col_pos - 1], seq2[row_pos - 1]))

            row_pos -= 1
            col_pos -= 1

    return (list(reversed(reverse_alignment)), alignment_matrix[-1, -1])
