"""Train the PAE aggregator to predict T cell antigen specificity from PAE matrices."""

import argparse
import logging
import os
import sys

import h5py
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from torch import nn, optim

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.models import PAEAggregator

logger = logging.getLogger()

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('training_data', help='path to the training data (HDF5 file)')
parser.add_argument('-o', '--output', required=True, help='path to output models to')
parser.add_argument('--seed', default=None, type=int, help='seed for random processes (Default: None)')
parser.add_argument('--label-name', default='label', help="name to use access label in dataset (Default: 'label')")
parser.add_argument('--num-epochs', type=int, default=5, help='number of epochs (Default: 5)')
parser.add_argument('--batch-size', type=int, default=32, help='training batch size (Default: 32)')
parser.add_argument('--learning-rate', type=float, default=1e-3, help='optimizer learning rate (Default: 1e-3)')
parser.add_argument(
    '--eval-batch-size', type=int, default=300, help='batch size for evaluating the model (Default: 300)'
)
parser.add_argument(
    '--eval-interval', type=int, default=50, help='Interval to evaluate model performance (Default: 50)'
)

add_logging_arguments(parser)


def evaluate_model(
    model: nn.Module,
    paes: np.ndarray,
    labels: np.ndarray,
    device: torch.device,
    batch_size: int = 1000,
) -> float:
    """Evaluate model using ROC-AUC.

    Sets model to evaluation mode.
    """
    model.eval()
    indices = np.arange(len(labels))

    num_eval_batches = int(np.ceil(len(indices) / batch_size))

    predictions = []
    for i in range(num_eval_batches):
        batch_idx_start = i * batch_size

        batch_indices = (
            indices[batch_idx_start : batch_idx_start + batch_size]
            if i < num_eval_batches - 1
            else indices[batch_idx_start:]
        )
        batch = torch.tensor(paes[batch_indices], dtype=torch.float32, device=device)

        prediction = model(batch)
        predictions.append(prediction.detach().cpu().numpy())

    predictions = np.concatenate(predictions)

    return roc_auc_score(labels, predictions)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level, args.log_file)

    if args.seed:
        logger.info('Setting random seed to %d', args.seed)
        torch.manual_seed(args.seed)
        random_generator = np.random.default_rng(args.seed)

    else:
        random_generator = np.random.default_rng()

    for argument, value in vars(args).items():
        logger.info('Parameter: %s=%r', argument, value)

    logger.info('Loading training data from %s', args.training_data)
    with h5py.File(args.training_data) as fh:
        paes = fh['pae'][:]
        labels = fh[args.label_name][:]
        folds = fh['fold'][:]

    indices = np.arange(len(labels))

    if not os.path.exists(args.output):
        logger.info('Creating output directory: %s', args.output)
        os.mkdir(args.output)

    for fold_idx in np.unique(folds):
        logger.info('Starting fold %d', fold_idx)

        training_indices = indices[folds != fold_idx]
        validation_indices = indices[folds == fold_idx]

        validation_data = paes[validation_indices]
        validation_labels = labels[validation_indices]

        logger.debug('Number of training data points: %d', len(training_indices))
        logger.debug('Number of validation data points: %d', len(validation_indices))

        logger.info('Initialising model')
        model = PAEAggregator(paes.shape[-1])
        model.train()

        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        logger.debug('Transferring model to %s', str(device))
        model.to(device)

        loss_fn = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=args.learning_rate)

        output_name = os.path.join(args.output, f'model_{fold_idx}.pt')

        batches_per_epoch = len(training_indices) // args.batch_size

        logger.info('Starting training')
        for epoch in range(1, args.num_epochs + 1):
            logger.info('Starting epoch %d', epoch)
            epoch_indices = training_indices.copy()
            random_generator.shuffle(epoch_indices)

            for i in range(batches_per_epoch):
                logger.debug('Training batch %d of %d', i + 1, batches_per_epoch)
                batch_indices = epoch_indices[i * args.batch_size : i * args.batch_size + args.batch_size]

                batch = torch.tensor(paes[batch_indices], dtype=torch.float32, device=device)
                batch_label = torch.tensor(labels[batch_indices], dtype=torch.float32, device=device)

                prediction = model(batch)

                loss = loss_fn(prediction, batch_label)
                loss.backward()

                optimizer.step()
                optimizer.zero_grad()

                if i % args.eval_interval == 0:
                    logger.info('Evaluating at iteration %d', (i + 1) + (batches_per_epoch * (epoch - 1)))

                    training_roc_auc = roc_auc_score(
                        batch_label.detach().cpu().numpy(), prediction.detach().cpu().numpy()
                    )

                    validation_roc_auc = evaluate_model(
                        model, validation_data, validation_labels, device, batch_size=args.eval_batch_size
                    )

                    logger.info('Loss: %f', loss.item())
                    logger.info('Training ROC-AUC: %f', training_roc_auc)
                    logger.info('Validation ROC-AUC: %f', validation_roc_auc)

                    model.train()

        logger.info('Saving model to %s', output_name)
        torch.save(model.state_dict(), output_name)


if __name__ == '__main__':
    main()
