#! /usr/bin/env python
"""Create data archive for the project."""

import json
import logging
import tarfile
import tempfile
from pathlib import Path

import numpy as np
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger

logger = logging.getLogger()
setup_logger(logger, level=snakemake.config['log_level'])


def add_confidences(tar: tarfile.TarFile, name: str, directory: Path, structures_name: str) -> None:
    """Add confidence scores to Tarball."""
    with open(str(directory / f'{name}_summary_confidences.json'), 'r') as fh:
        summary_confidence = json.load(fh)

    with open(str(directory / f'{name}_confidences.json'), 'r') as fh:
        confidence = json.load(fh)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir) / f'{name}_confidence.npz'
        np.savez_compressed(str(temp_path), **summary_confidence, **confidence, allow_pickle=False)

        tar.add(str(temp_path), arcname=f'predicted_structures/{structures_name}/{name}_confidence.npz')


logger.info('Creating output tarball at %s', snakemake.output[0])
with tarfile.open(snakemake.output[0], 'w:gz') as tar:
    logger.info('Adding sequences CSVs')
    tar.add(snakemake.input.sequences, arcname='sequences.csv')
    tar.add(snakemake.input.sequences_sample_50, arcname='sequences_sample_50.csv')
    tar.add(snakemake.input.sequences_sample_1500, arcname='sequences_sample_1500.csv')

    logger.info('Adding structures')
    tar.add(snakemake.input.structures_summary, arcname='structures/summary.csv')
    for path in Path(snakemake.input.structures).glob('*.pdb'):
        tar.add(str(path), arcname=f'structures/{path.name}')

    logger.info('Adding predicted structures')

    logger.info('Adding af3-af3-af3 structures')
    structures_name = 'model-af3/msa-af3/data-af3'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_af3_50):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_af3_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-mmseqs2-colabfold structures')
    structures_name = 'model-af3/msa-mmseqs2/data-colabfold'

    for path in (Path(p) for p in snakemake.input.structures_af3_mmseqs2_colabfold_50):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2 structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_50):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_paired_msas structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_paired_msas'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_paired_msas_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_paired_msas_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_unpaired_msas structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_unpaired_msas'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_templates structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_templates'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_templates_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_templates_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_unpaired_msas_no_paired_msas structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_unpaired_msas_no_paired_msas'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (
        Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_benchmark
    ):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_unpaired_msas_no_templates structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_unpaired_msas_no_templates'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_templates_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_templates_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_paired_msas_no_templates structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_paired_msas_no_templates'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_paired_msas_no_templates_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_paired_msas_no_templates_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-no_unpaired_msas_no_paired_msas_no_templates structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-no_unpaired_msas_no_paired_msas_no_templates'

    for path in (
        Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_no_templates_1500
    ):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (
        Path(p)
        for p in snakemake.input.structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_no_templates_benchmark
    ):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af3-af3-tcrmodel2-3_recycles structures')
    structures_name = 'model-af3/msa-af3/data-tcrmodel2/modifications-3_recycles'

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_3_recycles_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af3_af3_tcrmodel2_3_recycles_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding boltz-af3-tcrmodel2 structures')
    structures_name = 'model-boltz/msa-af3/data-tcrmodel2'

    for path in (Path(p) for p in snakemake.input.structures_boltz_af3_tcrmodel2_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_boltz_af3_tcrmodel2_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding chai-af3-tcrmodel2 structures')
    structures_name = 'model-chai/msa-af3/data-tcrmodel2'

    for path in (Path(p) for p in snakemake.input.structures_chai_af3_tcrmodel2_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_chai_af3_tcrmodel2_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_model.cif', ''), path.parent, structures_name)

    logger.info('Adding af2-af3-tcrmodel2 structures')
    structures_name = 'model-chai/msa-af3/data-tcrmodel2'

    for path in (Path(p) for p in snakemake.input.structures_af2_af3_tcrmodel2_1500):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_relaxed.pdb', ''), path.parent, structures_name)

    for path in (Path(p) for p in snakemake.input.structures_af2_af3_tcrmodel2_benchmark):
        tar.add(str(path), arcname=f'predicted_structures/{structures_name}/{path.name}')
        add_confidences(tar, path.name.replace('_relaxed.pdb', ''), path.parent, structures_name)

    logger.info('Adding PAEAggregator models')
    for path in (Path(p) for p in snakemake.input.PAEAggregator):
        tar.add(str(path), arcname=f'PAEAggregator/{path.name}')
