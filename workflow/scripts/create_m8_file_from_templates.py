#!/usr/bin/env python
"""Create m8 files from a templates directory."""

import glob
import hashlib
import json
import os
import re
from collections import OrderedDict, defaultdict

from Bio.PDB import MMCIFParser
from snakemake.script import snakemake


def create_cigar(query_indices: list[int], template_indices: list[int]) -> str:
    """Create MMSeqs2 style cigar representing the alignment from paired aligned residues."""
    query_idx = -1
    template_idx = -1

    cigar = []
    for next_query_idx, next_template_idx in zip(query_indices, template_indices, strict=True):
        query_step = next_query_idx - query_idx
        if query_step > 1:
            cigar.extend(['I' for _ in range(query_step - 1)])

        template_step = next_template_idx - template_idx
        if template_step > 1:
            cigar.extend(['D' for _ in range(template_step - 1)])

        cigar.append('M')
        query_idx = next_query_idx
        template_idx = next_template_idx

    return cigar


def format_cigar(cigar: list[str]) -> str:
    """Format cigar to condensed string representation (eg: 12I2D34M1I54M)."""
    count = 0
    prev_entry = None
    cigar_formatted = []
    for i, entry in enumerate(cigar):
        if entry == prev_entry:
            count += 1

        else:
            if i != 0:
                cigar_formatted.append((count, prev_entry))

            count = 1
            prev_entry = entry

    cigar_formatted.append((count, prev_entry))

    return ''.join(f'{count}{code}' for count, code in cigar_formatted)


with open(snakemake.input.msa_input, 'r') as fh:
    msa_input = json.load(fh)

sequence = msa_input['sequences'][0]['protein']['sequence']
sequence_id = hashlib.sha256(sequence.encode()).hexdigest()

template_paths = glob.glob(os.path.join(snakemake.input.templates, '*.cif'))
indices_paths = glob.glob(os.path.join(snakemake.input.templates, '*_indices.json'))

template_entries = defaultdict(dict)
for indice_path in indices_paths:
    entry_name = re.search(r'(\d+)_indices.json', os.path.basename(indice_path)).group(1)
    template_entries[entry_name]['indices_path'] = indice_path

for template_path in template_paths:
    entry_name = re.search(r'(\d+).cif', os.path.basename(template_path)).group(1)
    template_entries[entry_name]['structure_path'] = template_path

template_entries = sorted(template_entries.items())

m8_entries = []
for name, entry in template_entries:
    m8_entry = OrderedDict()

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('', entry['structure_path'])
    chain_name = next(structure[0].get_chains()).id

    with open(entry['indices_path'], 'r') as fh:
        indices = json.load(fh)

    cigar = create_cigar(indices['queryIndices'], indices['templateIndices'])

    m8_entry['query_id'] = sequence_id
    m8_entry['subject_id'] = f'{snakemake.wildcards.name}-{name}_{chain_name}'
    m8_entry['pident'] = ''
    m8_entry['length'] = ''
    m8_entry['mismatch'] = ''
    m8_entry['gapopen'] = ''
    m8_entry['query_start'] = str(min(indices['queryIndices']) + 1)
    m8_entry['query_end'] = str(max(indices['queryIndices']) + 1)
    m8_entry['subject_start'] = str(min(indices['templateIndices']) + 1)
    m8_entry['subject_end'] = str(max(indices['templateIndices']) + 1)
    m8_entry['evalue'] = ''
    m8_entry['bitscore'] = ''
    m8_entry['comment'] = format_cigar(cigar)

    m8_entries.append(m8_entry)

with open(snakemake.output[0], 'w') as fh:
    for entry in m8_entries:
        fh.write('\t'.join(entry.values()))
        fh.write('\n')
