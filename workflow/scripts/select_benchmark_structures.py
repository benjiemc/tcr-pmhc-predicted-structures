#!/usr/bin/env python
"""Script to select holdout structuers to evaluate AF3 models on."""

import os

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake

benchmark_set = pd.read_csv(snakemake.input.shi_parks_smith_benchmark)

benchmark_set.columns = [
    column.replace('(Å)', '').lower().replace(' ', '_').strip('_') for column in benchmark_set.columns
]
benchmark_set['pdb_id'] = benchmark_set['pdb_id'].str.lower()
benchmark_set['mhc_class'] = benchmark_set['mhc_class'].map({'Class I': 'MH1', 'Class II': 'MH2'})
benchmark_set['release_date'] = pd.to_datetime(benchmark_set['release_date'], format=r'%m/%d/%y')

benchmark_set = benchmark_set.rename({'pdb_id': 'pdb', 'mhc_class': 'mhc_type'}, axis='columns')

stcrdab_summary = pd.read_csv(os.path.join(snakemake.input.stcrdab, 'db_summary.dat'), delimiter='\t')
benchmark_set_annotated = benchmark_set.merge(
    stcrdab_summary[[col for col in stcrdab_summary.columns if col not in ('resolution', 'date')]],
    how='inner',
    on=['pdb', 'mhc_type'],
)


def is_constant(series: pd.Series) -> bool:
    """Return True if all values in the pandas Series are the same (or empty)."""
    array = series.to_numpy()
    return array.shape[0] == 0 or (array[0] == array).all()


def drop_constant_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a new DataFrame with constant columns removed."""
    return df.loc[:, ~df.apply(is_constant)]


benchmark_set_annotated = benchmark_set_annotated.dropna(axis='columns', how='all')
benchmark_set_annotated = drop_constant_columns(benchmark_set_annotated)

print(
    'Structures with annotation issues in STCRDab:',
    benchmark_set[~benchmark_set['pdb'].isin(benchmark_set_annotated['pdb'])]['pdb'].tolist(),
)

structure_summary = pd.read_csv(snakemake.input.structures_summary)
benchmark_set_annotated_processed = benchmark_set_annotated.merge(
    structure_summary,
    how='inner',
    on=['pdb', 'Achain', 'Bchain', 'antigen_chain', 'mhc_chain1', 'mhc_chain2', 'mhc_type'],
)
print(
    'Structure not selected for processing:',
    benchmark_set_annotated[~benchmark_set_annotated['pdb'].isin(benchmark_set_annotated_processed['pdb'])][
        'pdb'
    ].tolist(),
)


benchmark_set_annotated_processed = benchmark_set_annotated_processed.drop_duplicates('pdb', keep='first')


def get_sequences(path: str) -> dict[str, str]:
    """Get sequences for all chain in protein."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(
        os.path.basename(path).replace('.pdb', ''), os.path.join(snakemake.input.structures, path)
    )

    return {
        chain.id: ''.join(
            [IUPACData.protein_letters_3to1[res.get_resname().title()] for res in chain if res.id[0] == ' ']
        )
        for chain in structure[0]
    }


sequences = benchmark_set_annotated_processed['path'].apply(get_sequences)

filled_sequences = benchmark_set_annotated_processed.filter(regex=r'chain\d?(_x)?$').apply(
    lambda row: row.map(sequences[row.name].get), axis=1
)
filled_sequences = filled_sequences.dropna(axis='columns', how='all')
filled_sequences.columns = [col_name.replace('_x', '') + '_sequence' for col_name in filled_sequences.columns]

benchmark_set_annotated_processed = benchmark_set_annotated_processed.merge(
    filled_sequences, how='inner', right_index=True, left_index=True
)

benchmark_set_annotated_processed['id'] = benchmark_set_annotated_processed['path'].str.replace('.pdb', '')
benchmark_set_annotated_processed = benchmark_set_annotated_processed.set_index('id').drop('path', axis='columns')

benchmark_set_annotated_processed.to_csv(snakemake.output[0])
