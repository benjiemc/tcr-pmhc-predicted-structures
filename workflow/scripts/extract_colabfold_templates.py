#!/usr/bin/env python
"""Script for extracting template information from m8 encoded files."""

import gzip
import json
import os
import re
import tempfile

from Bio.PDB import MMCIFIO, Chain, MMCIFParser, Model, Polypeptide, Structure
from snakemake.script import snakemake

with open(snakemake.input[0]) as fh:
    templates = [line.rstrip().split() for line in fh if line]

cigars = [entry[12] for entry in templates]
ids = [entry[1].split('_') for entry in templates]

os.mkdir(snakemake.output[0])
parser = MMCIFParser(QUIET=True)
for i, ((pdb_id, chain_id), cigar) in enumerate(
    zip(ids[: snakemake.params.number_of_templates], cigars[: snakemake.params.number_of_templates], strict=True)
):
    # Extract Structural template
    path = os.path.join(snakemake.params.colabfold_data_dir, 'pdb', 'divided', pdb_id[1:3], pdb_id + '.cif.gz')

    with gzip.open(path, 'rt') as fh:
        structure = parser.get_structure(pdb_id, fh)
        date = parser._mmcif_dict['_pdbx_audit_revision_history.revision_date'][0]  # noqa: SLF001 workaround for dates

    new_structure = Structure.Structure(pdb_id + '_' + chain_id)
    new_model = Model.Model(0)
    new_chain = Chain.Chain('A')

    for residue in structure[0][chain_id]:
        if Polypeptide.is_aa(residue):
            new_chain.add(residue.copy())

    new_model.add(new_chain)
    new_structure.add(new_model)

    with tempfile.NamedTemporaryFile(suffix='.cif') as temp_file:
        io = MMCIFIO()
        io.set_structure(new_structure)
        io.save(temp_file.name)  # Will be overwritten but needed to produce the mmcif dictionary
        # Hack to set all chain IDs to 'A'
        mmcif_dict = io.dic
        mmcif_dict['_atom_site.label_asym_id'] = ['A' for _ in io.dic['_atom_site.label_asym_id']]

        # Hack to set valid seq_ids for all entries
        idx = 0
        new_seq_ids = []
        prev_seq_id = None
        for seq_id in mmcif_dict['_atom_site.label_seq_id']:
            if seq_id != prev_seq_id:
                idx += 1
                prev_seq_id = seq_id

            new_seq_ids.append(str(idx))

        mmcif_dict['_atom_site.label_seq_id'] = new_seq_ids

        io.set_dict(mmcif_dict)
        io.save(temp_file.name)

        with open(temp_file.name, 'r') as fh:
            lines = fh.readlines()

    # Hack to add the date fields required by AF3
    lines.insert(1, f'#\n_pdbx_audit_revision_history.revision_date {date}\n')

    with open(os.path.join(snakemake.output[0], f'{i}.cif'), 'w') as fh:
        fh.write(''.join(lines))

    # Get query mapping
    indices = {
        'queryIndices': [],
        'templateIndices': [],
    }

    query_idx = 0
    template_idx = 0
    for repeats, action in re.findall(r'(?P<number>\d+)(?P<action>[IMD]+)', cigar):
        for _ in range(int(repeats)):
            match action:
                case 'M':
                    indices['queryIndices'].append(query_idx)
                    indices['templateIndices'].append(template_idx)

                    query_idx += 1
                    template_idx += 1

                case 'I':
                    query_idx += 1

                case 'D':
                    template_idx += 1

                case _:
                    msg = f'Invalid action: {action}'
                    raise ValueError(msg)

    with open(os.path.join(snakemake.output[0], f'{i}_indices.json'), 'w') as fh:
        json.dump(indices, fh)
