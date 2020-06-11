"""The CAT logger.

Index
-----
.. currentmodule:: CAT.logger
.. autosummary::
    logger

API
---
.. autodata:: logger
    :annotation: = logging.RootLogger

"""

import logging

from scm.plams import Molecule
from scm.plams.core.basejob import Job

__all__ = ['logger']

#: A logger for CAT
logger = logging.getLogger('CAT')
logging.basicConfig(level=logging.DEBUG,
                    format='[%(asctime)s] %(levelname)s: %(message)s',
                    datefmt='%H:%M:%S')


def log_start(job: Job, mol: Molecule, job_preset: str, job_name: str) -> None:
    """Log the start of a PLAMS :class:`Job`."""
    msg = 'has started'
    args = job.__class__.__name__, mol.properties.name, job_preset, job_name, msg
    logger.info('{}: {} {} ({}) {}'.format(*args))


def log_succes(job: Job, mol: Molecule, job_preset: str, job_name: str) -> None:
    """Log the succesful termination of a PLAMS :class:`Job`."""
    msg = 'is successful'
    args = job.__class__.__name__, mol.properties.name, job_preset, job_name, msg
    logger.info('{}: {} {} ({}) {}'.format(*args))


def log_fail(job: Job, mol: Molecule, job_preset: str, job_name: str) -> None:
    """Log the unsuccesful termination of a PLAMS :class:`Job`."""
    msg = 'has failed'
    args = job.__class__.__name__, mol.properties.name, job_preset, job_name, msg
    logger.info('{}: {} {} ({}) {}'.format(*args))


def log_copy(job: Job, mol: Molecule, job_preset: str, job_name: str, job_old: Job) -> None:
    """Log the copying of a previous PLAMS :class:`Job`."""
    msg = f'results copied from {job_old.name}'
    args = job.__class__.__name__, mol.properties.name, job_preset, job_name, msg
    logger.info('{}: {} {} ({}) {}'.format(*args))
