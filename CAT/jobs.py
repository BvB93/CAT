"""A module designed for running Jobs.

Index
-----
.. currentmodule:: CAT.jobs
.. autosummary::
    get_main_molecule
    _xyz_to_mol
    get_energy
    _get_name
    pre_process_settings
    retrieve_results
    job_single_point
    job_geometry_opt
    job_freq

API
---
.. automethod:: get_main_molecule
.. autofunction:: _xyz_to_mol
.. automethod:: get_energy
.. autofunction:: _get_name
.. autofunction:: pre_process_settings
.. autofunction:: retrieve_results
.. autofunction:: job_single_point
.. autofunction:: job_geometry_opt
.. autofunction:: job_freq

"""

from shutil import rmtree
from typing import (Optional, Type)
from os.path import join

import numpy as np

from scm.plams.core.basejob import Job
from scm.plams import (Molecule, Settings, Results, config, add_to_class, ResultsError,
                       ADFJob, AMSJob, Units, Cp2kResults)

import qmflows

from .logger import (logger, log_start, log_succes, log_fail, log_copy)
from .thermo_chem import get_thermo
from .utils import type_to_string

__all__ = ['job_single_point', 'job_geometry_opt', 'job_freq']


@add_to_class(Cp2kResults)
def get_main_molecule(self) -> Optional[Molecule]:
    try:
        filename = self.files['cp2k-pos-1.xyz']
        return _xyz_to_mol(filename)
    except TypeError:
        pass

    iterator = (join(self.job.path, f) for f in self.files if '.xyz' in f)
    for filename in iterator:
        return _xyz_to_mol(filename)
    raise ResultsError(f'Failed to retrieve main molecule from {self.job.name}')


def _xyz_to_mol(filename: str) -> Molecule:
    """Grab the last geometry from an .xyz file and return it as a :class:`Molecule` instance."""
    with open(filename, 'r') as f:
        atom_count = int(next(f))
        for line_count, _ in enumerate(f, 2):
            pass

    mol_count = line_count // (2 + atom_count)
    return Molecule(filename, geometry=mol_count)


@add_to_class(Cp2kResults)
def get_energy(self, index: int = -1, unit: str = 'Hartree') -> float:
    """Return the energy of the last occurence of ``'ENERGY| Total FORCE_EVAL'`` in the output."""
    energy_str = self.grep_output('ENERGY| Total FORCE_EVAL')[index]
    energy = float(energy_str.rsplit(maxsplit=1)[1])
    return Units.convert(energy, 'Hartree', unit)


def _get_name(name: str) -> str:
    manager = config.default_jobmanager
    if name in manager.names:
        number = str(1 + manager.names[name]).zfill(manager.settings.counter_len)
        return f'{name}.{number}'
    return name


def pre_process_settings(mol: Molecule, s: Settings,
                         job_type: Type[Job], template_name: str) -> Settings:
    """Update all :class:`Settings`, **s**, with those from a QMFlows template (see **job**)."""
    ret = Settings()
    type_key = type_to_string(job_type)
    ret.input = getattr(qmflows, template_name)['specific'][type_key].copy()
    ret.update(s)

    if job_type is AMSJob:
        mol.properties.pop('charge', None)
        # ret.input.ams.system.BondOrders._1 = adf_connectivity(mol)
        if 'uff' not in s.input:
            ret.input.ams.system.charge = int(sum(
                at.properties.get('charge', 0) for at in mol
            ))

    elif job_type is ADFJob:
        mol.properties.pop('charge', None)
        if not ret.input.charge:
            ret.input.charge = int(sum(at.properties.get('charge', 0) for at in mol))
    return ret


JOB_PRESETS = frozenset({
    'geometry optimization', 'frequency analysis', 'single point', 'MD calculation'
})


def retrieve_results(mol: Molecule, results: Results, job_preset: str) -> None:
    """Unpack the :class:`results` from a PLAMS-facilitated calculation.

    Performs an inplace update of **results**:

        * The Cartesian coordinates of :attr:`Results.job.molecule` are updated.
        * Energies and frequencies are added to :attr:`Results.job.molecule.properties`.
        *

    Paramaters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    results : |plams.Results|_
        A PLAMS :class:`Results` instance.

    job_preset : str
        The name of a job preset from :mod:`CAT.jobs`.
        Accepted values are ``"single point"``, ``"geometry optimization"`` or
        ``"frequency analysis"``.

    """
    if job_preset not in JOB_PRESETS:
        raise ValueError(f'Invalid value for job_preset: {job_preset!r}')

    # Unpack arguments
    job = results.job
    name = job.name

    # Define more arguments
    freq = np.empty(1)
    nan_dict = {
        'frequencies': np.nan,
        'energy': {'E': np.nan, 'H': np.nan, 'S': np.nan, 'G': np.nan}
    }

    try:
        if job.status in {'failed', 'crashed'}:
            raise _get_results_error(results)

        # Read all relevant results
        energy = mol.properties.energy.E = results.get_energy(unit='kcal/mol')
        if job_preset in ('geometry optimization', 'frequency analysis'):
            mol_new = results.get_main_molecule() or Molecule()
            mol.from_array(mol_new)

        if job_preset == 'frequency analysis':
            freq = mol.properties.frequencies = results.get_frequencies()
            energy = mol.properties.energy = get_thermo(mol, freq, energy)

        # Evaluate all results
        if not (energy and isinstance(freq, np.ndarray)):
            raise _get_results_error(results)
        log_succes(job, mol, job_preset, name)

    except Exception as ex:  # Failed to retrieve results
        if job_preset == 'geometry optimization':
            mol.properties.is_opt = False
        mol.properties.soft_update(nan_dict)
        log_fail(job, mol, job_preset, name)
        logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)

    else:
        if job_preset == 'geometry optimization':
            mol.properties.is_opt = True
        if job.status != 'copied':
            return None

        # results.job is a copy from a previously run job.
        # Exchange the attributes of results.job wth those of the revious job
        job_old = config.default_jobmanager.hashes[job._hash]()
        log_copy(job, mol, job_preset, name, job_old)
        rmtree(job.path)
        for key, value in vars(job_old).items():
            if key != 'molecule':
                setattr(job, key, value)
    return None


def _get_results_error(results: Results) -> ResultsError:
    """Raise a :exc:`ResultsError` with the content of ``results['$JN.err']`` as error mesage."""
    filename = results['$JN.err']
    with open(filename, 'r') as f:
        iterator = (i.rstrip('\n') for i in f)
        for item in iterator:
            if item:
                return ResultsError(item)
        else:
            return ResultsError()


@add_to_class(Job)
def _finalize(self):
    """Modified PLAMS :meth:`Job._finalize` method.

    1. All references to the :func:`scm.plams.log` function are removed.
    2. A $JN.hash file is created which contains the hash of this instance.

    Gather the results of the job execution and organize them.
    This method collects steps 9-12 from :ref:`job-life-cycle`.
    Should not be overridden.

    """
    if config.preview is False:
        self.results.collect()
        self.results.finished.set()
        if self.status != 'crashed':
            self.status = 'finished'
            self._log_status(3)
            if self.check():
                self.results._clean(self.settings.keep)
                self.postrun()
                self.status = 'successful'
                if self.settings.pickle:
                    self.pickle()
                    filename = join(self.path, f'{self.name}.hash')
                    with open(filename, 'w') as f:
                        try:  # Will raise a TypeError if rerun prevetion is disabled
                            f.write(self.hash())
                        except TypeError:
                            pass
            else:
                self.status = 'failed'
    else:
        self.status = 'preview'
        self.results.finished.set()
    self.results.done.set()

    if self.parent and self in self.parent:
        self.parent._notify()

    self._log_status(1)


@add_to_class(Molecule)
def job_single_point(self, job_type: Type[Job],
                     settings: Settings,
                     name: str = 'Single_point',
                     ret_results: bool = False,
                     read_template: bool = True) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies.

    Paramaters
    ----------
    job_type : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    read_template : bool
        Whether or not to update **settings** using a QMFlows template or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Grab the default settings for a specific job and update them with user provided settings
    _sp = 'singlepoint'
    s = settings if not read_template else pre_process_settings(self, settings, job_type, _sp)

    # Run the job; extract energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'single point', _name)

    results = job.run()
    retrieve_results(self, results, 'single point')

    inp_name = join(job.path, f'{job.name}.in')
    try:
        self.properties.job_path.append(inp_name)
    except TypeError:
        self.properties.job_path = [inp_name]

    # Return results
    if ret_results:
        return job.results
    return None


@add_to_class(Molecule)
def job_geometry_opt(self, job_type: Type[Job],
                     settings: Settings,
                     name: str = 'Geometry_optimization',
                     ret_results: bool = False,
                     read_template: bool = True) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies and final geometries.

    Paramaters
    ----------
    job_type : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    read_template : bool
        Whether or not to update **settings** using a QMFlows template or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Grab the default settings for a specific job and update them with user provided settings
    _geo = 'geometry'
    s = settings if not read_template else pre_process_settings(self, settings, job_type, _geo)

    # Run the job; extract geometries and energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'geometry optimization', _name)

    results = job.run()
    retrieve_results(self, results, 'geometry optimization')

    inp_name = join(job.path, f'{job.name}.in')
    try:
        self.properties.job_path.append(inp_name)
    except TypeError:
        self.properties.job_path = [inp_name]

    # Return results
    if ret_results:
        return job.results
    return None


@add_to_class(Molecule)
def job_md(self, job_type: Type[Job],
           settings: Settings,
           name: str = 'MD',
           ret_results: bool = False,
           read_template: bool = False) -> Optional[Results]:
    """Function for running an arbritrary jobs, extracting total energies and final geometries.

    Paramaters
    ----------
    job_type : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    read_template : bool
        Whether or not to update **settings** using a QMFlows template or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Grab the default settings for a specific job and update them with user provided settings
    s = settings

    # Run the job; extract geometries and energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'MD calculation', _name)

    results = job.run()
    retrieve_results(self, results, 'MD calculation')

    inp_name = join(job.path, f'{job.name}.in')
    try:
        self.properties.job_path.append(inp_name)
    except TypeError:
        self.properties.job_path = [inp_name]

    # Return results
    if ret_results:
        return job.results
    return None


@add_to_class(Molecule)
def job_freq(self, job_type: Type[Job],
             settings: Settings,
             name: str = 'Frequency_analysis',
             opt: bool = True,
             ret_results: bool = False,
             read_template: bool = True) -> Optional[Results]:
    """Function for running an arbritrary Jobs.

    Extracts total energies, final geometries and
    thermochemical quantities derived from vibrational frequencies.

    Paramaters
    ----------
    job : |Callable|_
        A type Callable of a class derived from :class:`Job`, e.g. :class:`AMSJob`
        or :class:`Cp2kJob`.

    settings : |plams.Settings|_
        The settings for **job**.

    name : str
        The name of **job**.

    opt : bool
        Perform a geometry optimization (see :func:`.job_geometry_opt`) before calculating
        frequencies.

    ret_results : bool
        Whether or not the :class:`Results` instance should be returned or not.

    read_template : bool
        Whether or not to update **settings** using a QMFlows template or not.

    Returns
    -------
    |plams.Results|_
        Optional: If ``ret_results=True` return the :class:`Results` instance produced by this job.

    """
    # Preceed the frequency analysis with a geometry optimization
    if opt:
        self.job_geometry_opt(job_type, settings)

    # Grab the default settings for a specific job and update them with user provided settings
    s = settings if not read_template else pre_process_settings(self, settings, job_type, 'freq')

    # Run the job; extract geometries and (Gibbs free) energies
    job = job_type(molecule=self, settings=s, name=name)
    _name = _get_name(name)
    log_start(job, self, 'frequency analysis', _name)

    results = job.run()
    retrieve_results(self, results, 'frequency analysis')

    inp_name = join(job.path, f'{_name}.in')
    try:
        self.properties.job_path.append(inp_name)
    except TypeError:
        self.properties.job_path = [inp_name]

    # Return results
    if ret_results:
        return job.results
    return None
