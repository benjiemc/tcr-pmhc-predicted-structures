import logging
from argparse import ArgumentParser

LOG_MAP = {
    'error': logging.ERROR,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG,
}


def add_logging_arguments(parser: ArgumentParser) -> None:
    """Add logging arguments to parser."""
    logging_group = parser.add_argument_group('Logging', 'Options for logging')
    logging_group.add_argument(
        '--log-level',
        choices=['debug', 'info', 'warning', 'error'],
        default='warning',
        help="Level to log messages at (Default: 'warning')",
    )
    logging_group.add_argument('--log-file', default=None, help='File to output logs (default is to write to stderr)')


def setup_logger(logger, level='warning'):
    """Set up logger for command line applications.

    Args:
        logger: logger object from logging.getLogger(__name__) call
        level: one of 'error', 'warning', 'info', or 'debug'

    """
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logger.addHandler(handler)

    logger.setLevel(LOG_MAP[level])
