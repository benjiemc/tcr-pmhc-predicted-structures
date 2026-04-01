"""Run PAEAggregator model on predictions."""

import argparse
import glob
import logging
import os
import sys

import h5py
import numpy as np
import pandas as pd
import torch

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.models import PAEAggregator

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    '--mode', choices=['untrained', 'evaluation', 'run'], default='run', help="prediction mode (Default: 'run')"
)

inputs = parser.add_argument_group('Inputs')
inputs.add_argument('--paes', required=True, help='path to pae dataset')
inputs.add_argument('--model', help='path to model')

outputs = parser.add_argument_group('Outputs')
outputs.add_argument('--output', '-o', required=True, help='path to output predictions in csv format')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading PAEs from %s', args.paes)
    with h5py.File(args.paes, 'r') as fh:
        paes = fh['pae'][:]

        names = np.array([name.decode() for name in fh['id'][:]])
        seeds = fh['seed'][:]
        samples = fh['sample'][:]

        if args.mode == 'evaluation':
            folds = fh['fold'][:]

    model_size = paes.shape[1]
    logger.debug('Model size is %d', model_size)

    logger.info('Making predictions')
    match args.mode:
        case 'untrained':
            model = PAEAggregator(model_size)
            model.eval()

            predictions = model(torch.Tensor(paes)).cpu().detach().numpy()

            predictions_df = pd.DataFrame(
                {
                    'name': names,
                    'seed': seeds,
                    'sample': samples,
                    'prediction': predictions,
                }
            )

        case 'evaluation':
            predictions_df = pd.DataFrame()
            for fold in sorted(np.unique(folds)):
                model_path = os.path.join(args.model, f'model_{fold}.pt')
                model = PAEAggregator(model_size)

                logger.debug('Loading model from %s', model_path)
                model.load_state_dict(torch.load(model_path, weights_only=True, map_location=torch.device('cpu')))
                model.eval()

                fold_mask = folds == fold

                predictions = model(torch.Tensor(paes[fold_mask])).cpu().detach().numpy()

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

        case 'run':
            predictions_df = pd.DataFrame({'name': names, 'seed': seeds, 'sample': samples})
            for model_path in sorted(glob.glob(os.path.join(args.model, 'model_*.pt'))):
                model_name = os.path.basename(model_path).replace('model_', '').replace('.pt', '')

                model = PAEAggregator(model_size)

                logger.debug('Loading model from %s', model_path)
                model.load_state_dict(torch.load(model_path, weights_only=True, map_location=torch.device('cpu')))
                model.eval()

                predictions = model(torch.Tensor(paes)).cpu().detach().numpy()
                predictions_df[f'prediction_{model_name}'] = predictions

    logger.info('Outputting predictions to %s', args.output)
    predictions_df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
