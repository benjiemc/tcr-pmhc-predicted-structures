#!/usr/bin/env python
"""Run PAEAggregator model on predictions."""

import logging
import os

import h5py
import numpy as np
import pandas as pd
import torch
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.models import PAEAggregator

logger = logging.getLogger()
setup_logger(logger, snakemake.config['log_level'])

logger.info('Loading PAEs from %s', snakemake.input.paes)
with h5py.File(snakemake.input.paes, 'r') as fh:
    paes = fh['pae'][:]

    names = np.array([name.decode() for name in fh['id'][:]])
    seeds = fh['seed'][:]
    samples = fh['sample'][:]

    folds = fh['fold'][:]

logger.info('Loading scaling data from %s', snakemake.input.pcas)
scaling_data_df = pd.read_csv(snakemake.input.pcas)

scaling_data_df['name'] = scaling_data_df['name'].astype(str)
scaling_data_df = scaling_data_df.set_index(['name', 'seed', 'sample'])

scaling_data = np.stack(
    [
        scaling_data_df.loc[name, seed, sample].to_numpy()
        for name, seed, sample in zip(names, seeds, samples, strict=True)
    ]
)

logger.info('Making predictions')
match snakemake.wildcards.training_status:
    case 'untrained':
        model = PAEAggregator(snakemake.params.model_size)
        model.eval()

        predictions = model(torch.Tensor(paes), torch.Tensor(scaling_data)).cpu().detach().numpy()
        latents = torch.mul(model.weights, torch.Tensor(paes)).cpu().detach().numpy()

        predictions_df = pd.DataFrame(
            {
                'name': names,
                'seed': seeds,
                'sample': samples,
                'prediction': predictions,
            }
        )

    case 'trained':
        predictions_df = pd.DataFrame()
        for fold in sorted(np.unique(folds)):
            model_path = os.path.join(snakemake.input.model, f'model_{fold}.pt')
            model = PAEAggregator(snakemake.params.model_size)

            logger.debug('Loading model from %s', model_path)
            model.load_state_dict(torch.load(model_path, weights_only=True, map_location=torch.device('cpu')))
            model.eval()

            fold_mask = folds == fold

            predictions = (
                model(torch.Tensor(paes[fold_mask]), torch.Tensor(scaling_data[fold_mask])).cpu().detach().numpy()
            )

            predictions_df = pd.concat(
                [
                    predictions_df,
                    pd.DataFrame(
                        {
                            'name': names[fold_mask],
                            'seed': seeds[fold_mask],
                            'sample': samples[fold_mask],
                            'prediction': predictions,
                        }
                    ),
                ]
            )

logger.info('Outputting predictions to %s', snakemake.output[0])
predictions_df.to_csv(snakemake.output[0], index=False)
