"""A module to hold the :class:`GenJobManager` class.

Index
-----
.. currentmodule:: CAT.gen_job_manager
.. autosummary::
    GenJobManager

API
---
.. autoclass:: GenJobManager
   :members:
   :private-members:
   :special-members:

"""

import os
from typing import (Optional, Callable)
from os.path import (join, isfile, abspath, isdir, exists, normpath)

import dill as pickle

from scm.plams import (JobManager, FileError, Settings, PlamsError)
from scm.plams.core.basejob import Job

__all__ = ['GenJobManager']


class GenJobManager(JobManager):
    """A modified version of the PLAMS :class:`JobManager` class.

    In order to minimize its memory footprint, :meth:`GenJobManager.load_job` will now populate
    the values of :attr:`GenJobManager.hashes` with a callable the creates and returns a
    :class:`Job` instance (rather than permantly sotring the :class:`Job` instance).

    Attributes
    ----------
    foldername : str
        The name of the working folder (see :attr:`GenJobManager.workdir`).

    workdir : str
        The absolute path to the working folder.

    logfile : str
        The absolute path+filename of the logfile.

    input : str
        The absolute path to the copy of the input file in the working folder.

    settings : |plams.Settings|_
        A :class:`Settings` instance with settings for this job manager.

    jobs : |list|_ [|plams.Jobs|_]
        A list of all jobs managed with this instance (in order of run calls).

    names : |dict|_ [|str|_, |int|_]
        A dictionary with names of jobs.
        For each name an integer value is stored indicating
        how many jobs with that basename have already been run.

    hashes : |dict|_ [|str|_, |Callable|_]
        A dictionary working as a hash-table for jobs.
        Values created by :meth:`GenJobManager.load_job` are stored as callables which
        in turn create :class:`Job` instances.

    """

    def __init__(self, settings: Settings,
                 path: Optional[str] = None,
                 folder: Optional[str] = None,
                 hashing: str = 'input') -> None:
        """Initialize the :class:`GenJobManager` instance."""
        self.settings = settings
        self.jobs = []
        self.names = {}
        self.hashes = {}

        if path is None:
            self.path = os.getcwd()
        elif isdir(path):
            self.path = abspath(path)
        else:
            raise PlamsError(f'Invalid path: {path}')

        basename = normpath(folder) if folder else 'plams_workdir'
        self.foldername = basename
        self.workdir = join(self.path, self.foldername)
        self.logfile = join(self.workdir, 'logfile')
        self.input = join(self.workdir, hashing)
        if not exists(self.workdir):
            os.mkdir(self.workdir)

    @staticmethod
    def _unpickle(filename: str) -> Optional[Job]:
        """Attempt to unpickle and return a :class:`Job` containing file."""
        with open(filename, 'rb') as f:
            try:
                return pickle.load(f)
            except Exception:
                return None

    def _get_job(self, filename: str) -> Callable:
        """Return a callable which converts **filename** into a :class:`Job` instance."""
        def unpickle_job() -> Optional[Job]:
            _filename = filename.replace('.hash', '.dill')
            ret = GenJobManager._unpickle(_filename)
            ret.jobmanager = self
            return ret
        return unpickle_job

    def load_job(self, filename: str) -> None:
        """Load a previously saved job from **filename**, populating :attr:`GenJobManager.hashes`.

        The hash of the loaded job is stored as key in :attr:`GenJobManager.hashes`,
        the matching value being a callable that loads the actual the :class:`Job` instance.

        Parameters
        ----------
        filename : str
            A path to a .hash file in some job folder.
            A :class:`Job` instance stored there is loaded and returned.
            All attributes of this instance removed before pickling are restored.

        """
        # Raise an error if **filename** cannot be found
        if isfile(filename):
            filename = abspath(filename)
        else:
            raise FileError(f'File {filename} not present')

        # Read the has from filename
        with open(filename, 'r') as f:
            h = f.read().rstrip('\n')
        self.hashes[h] = self._get_job(filename)

    def remove_job(self, job: Job) -> None:
        """Remove **job** from the job manager; forget its hash."""
        if job in self.jobs:
            self.jobs.remove(job)
            job.jobmanager = None
        h = job.hash()
        if h in self.hashes:
            del self.hashes[h]

    def _check_hash(self, job: Job) -> Optional[Job]:
        """Calculate and check the hash of **job**.

        If the hash is not ``None``, previous jobs are searched for the same hash.
        Matches are returned if applicable; returns ``None`` otherwise.

        Parameter
        ---------
        job : |plams.Job|_
            A :class:`Job` instance.

        Returns
        -------
        |plams.Job|_
            Optional: A matching :class:`Job` instance stored in :attr:`GenJobManager.hashes`.
            If no matching hash is found return ``None``.

        """
        h = job.hash()
        setattr(job, '_hash', h)
        if h is None:  # No hash available for **job**, move along
            return None

        if h in self.hashes:
            try:
                return self.hashes[h]()
            except AttributeError:  # In case the job unpickling fails
                pass

        filename = join(job.path, job.name)
        func = self._get_job(filename + '.dill')
        self.hashes[h] = func  # Set a callable that returns a Job instance
        return None
