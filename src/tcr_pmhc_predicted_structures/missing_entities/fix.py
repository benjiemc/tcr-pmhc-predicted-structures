"""Functions and classes for modelling missing residues in structures using MODELLER.

Requires modeller to be available can be downloaded here: https://salilab.org/modeller/

"""

import contextlib
import io
import logging
import os
import sys
from tempfile import TemporaryDirectory

from Bio.PDB import PDBIO, Chain, PDBParser, Structure

logger = logging.getLogger(__name__)

try:
    from modeller import Environ, Selection, log
    from modeller.automodel import AutoModel

except ImportError:
    logger.exception('Modeller not found. Please follow install instructions here: https://salilab.org/modeller/')
    sys.exit(1)

_ALIGNMENT_TEMPLATE = '''>P1;target
sequence:target:::::::0.00:0.00
{sequence}*

>P1;struct
structure:struct:{start}:{chain_id}:::::0.00:0.00
{template}*'''


class MissingResiduesModel(AutoModel):
    """AutoModel class for modelling missing resiudes.

    Only the missing residues can be optimized.

    Attributes:
        missing_residue_ranges: list of ranges (inclusive) of residues missing in the template.

    """

    def __init__(
        self,
        missing_residue_selection: list[str],
        missing_atom_selection: list[str],
        *args: tuple,
        **kwargs: dict,
    ) -> None:
        """Initialize the class.

        Args:
            missing_residue_selection: paths to the missing residues in MODELLER format (eg: '32:A')
            missing_atom_selection: paths to the missing atoms in MODELLER format (eg: 'CD1:24:A')
            *args: arguments for AutoModel class
            **kwargs: key-worded arguments for AutoModel class

        """
        super().__init__(*args, **kwargs)
        self.missing_residue_selection = missing_residue_selection
        self.missing_atom_selection = missing_atom_selection

    def select_atoms(self) -> Selection:
        """Atoms selected for optimization are only in the missing residues range."""
        selection = Selection()

        for missing_residue in self.missing_residue_selection:
            # Changing chain to A as required by MODELLER
            missing_residue_info = missing_residue.split(':')
            missing_residue_path = ':'.join([*missing_residue_info[:-1], 'A'])

            selection.add(self.residues[missing_residue_path])

        for missing_atom in self.missing_atom_selection:
            # Changing chain to A as required by MODELLER
            missing_atom_info = missing_atom.split(':')
            missing_atom_path = ':'.join([*missing_atom_info[:-1], 'A'])

            selection.add(self.atoms[missing_atom_path])

        return selection


def predict_missing_residues(
    alignment: list[tuple[str, str]],
    structure: Structure.Structure,
    missing_residue_selection: list[str],
    missing_atom_selection: list[str],
    chain_id: str,
) -> Chain.Chain:
    """Use MODELLER to fill in missing residues on a protein chain.

    Args:
        alignment: aligned complete sequence (including missing residues) and template sequence from structure
        structure: template PDB structure
        missing_residue_selection: paths to the missing residues in MODELLER format (eg: '32:A')
        missing_atom_selection: paths to the missing atoms in MODELLER format (eg: 'CD1:24:A')
        chain_id: ID of the chain to fill in missing residues

    Returns:
        new chain with filled in residues

    """
    sequence = ''.join([res for res, _ in alignment])
    template = ''.join([res for _, res in alignment])

    residues = [res for res in structure[0][chain_id].get_residues() if res.id[0] == ' ']
    start_pos = f'{residues[0].id[1]}{residues[0].id[2].strip()}'

    with TemporaryDirectory() as temp_dir:
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(temp_dir, 'struct.pdb'))

        align_path = os.path.join(temp_dir, 'alignment.ali')
        with open(align_path, 'w') as fh:
            fh.write(
                _ALIGNMENT_TEMPLATE.format(sequence=sequence, template=template, chain_id=chain_id, start=start_pos)
            )

        with contextlib.chdir(temp_dir):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                log.none()
                env = Environ()
                env.io.atom_files_directory = ['.']

                prediction = MissingResiduesModel(
                    missing_residue_selection,
                    missing_atom_selection,
                    env,
                    alnfile='alignment.ali',
                    knowns='struct',
                    sequence='target',
                    root_name='output_struct',
                )

                prediction.starting_model = 1
                prediction.ending_model = 1

                prediction.make()

            logger.debug(output.getvalue())

        pdb_parser = PDBParser(QUIET=True)
        fixed_structure = pdb_parser.get_structure('', os.path.join(temp_dir, 'output_struct.B99990001.pdb'))

    fixed_chain = fixed_structure[0]['A']
    fixed_chain.detach_parent()

    fixed_chain.id = chain_id

    return fixed_chain
