import os
import operator
from shutil import rmtree
from pathlib import Path
from contextlib import AbstractContextManager
from typing import (
    Optional, Union, Dict, Hashable, MutableMapping, TypeVar, Iterable, Container, Tuple, Callable,
    Any, List, Type
)

import numpy as np
import pandas as pd

import rdkit
import qmflows
from rdkit.Chem.AllChem import UFFGetMoleculeForceField as UFF
from scm.plams import finish, Settings, Molecule
from scm.plams.core.basejob import Job
from assertionlib.dataclass import AbstractDataClass

from ..utils import restart_init, _job_dict
from ..logger import logger
from ..settings_dataframe import SettingsDataFrame
from ..workflows.workflow_dicts import finalize_templates as load_templates
from ..frozen_settings import FrozenSettings

T = TypeVar('T')

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
MOL = ('mol', '')
OPT = ('opt', '')


def _return_True(value: Any) -> bool:
    """Return ``True``."""
    return True


def _lt_0(value) -> int:
    """Return if **value** is smaller than ``0``."""
    return value < 0


def pop_and_concatenate(mapping: MutableMapping[Hashable, T], base_key: Hashable,
                        filter_func: Callable[[Any], bool] = _return_True) -> Tuple[T, ...]:
    """Take a key and :meth:`pop<dict.pop>` all values from **mapping**.

    The popping will continue as long as :code:`base_key + str(i)` is available in the mapping,
    where ``i`` is an :class:`int` larger than 1.
    The value if ``i`` will start from 1 and increase by `+ 1` every iteration.

    Examples
    --------
    .. code:: python

        >>> mapping: dict = {
        ...     'job1': 1,
        ...     'job2': 2,
        ...     'job3': 3,
        ...     'final_key': True
        ... }

        >>> base_key: str = 'job'
        >>> value_tuple: tuple = concatenate_values(mapping, base_key)
        >>> print(value_tuple)
        (1, 2, 3)

        >>> print(mapping)
        {'final_key': True}

    Parameters
    ----------
    mapping : :data:`MutableMapping`
        A dictionary or other mutable mapping.

    base_key : :data:`Hashable<typing.Hashable>`
        The base key which will be appended with successively increasing integers.

    filter_func : :data:`Callable<typing.Callable>`
        A callable for truth-testing each extracted **mapping** value.
        Values returning ``False`` will not be added to the to-be returned :class:`tuple`.

    Returns
    -------
    :class:`tuple`
        A tuple with values popped from **mapping**.

    """
    i = 1
    ret = []
    while True:
        key = f'{base_key}{i}'
        try:
            value = mapping.pop(key)
        except KeyError:
            return tuple(ret)
        else:
            if filter_func(value):
                ret.append(value)
            i += 1


OptionalJobType = Union[None, Type[Job], Iterable[None], Iterable[Type[Job]]]
OptionalSettings = Union[None, Settings, Iterable[None], Iterable[Settings]]


class WorkFlow(AbstractDataClass):
    """A workflow manager.

    Examples
    --------
    Typical usage example:

    .. code:: python

        >>> import pandas as pd

        >>> # Prepare workflow parameters
        >>> df = pd.DataFrame(...)
        >>> settings = Settings(...)
        >>> def fancy_df_func(df, **kwargs):
        ...     pass

        >>> # Create the workflow
        >>> workflow = WorkFlow.from_template(settings, name='asa')
        >>> print(workflow)
        WorkFlow(
            name       = 'asa',
            db         = None,
            read       = False,
            write      = False,
            overwrite  = False,
            path       = '.',
            keep_files = True,
            jobs       = None,
            settings   = None
        )

        >>> # Run the workflow
        >>> idx = workflow.from_db(df)
        >>> workflow(fancy_df_func, df, index=idx)

        >>> # Export all workflow results
        >>> job_recipe = workflow.get_recipe()
        >>> workflow.to_db(df, job_recip=job_recipe)

    """

    #: Map a name to a workflow template.
    _WORKFLOW_TEMPLATES: FrozenSettings = load_templates()

    #: A context manager for supressing Pandas :exc:`SettingwithCopyWarning`.
    _SUPRESS_SETTINGWITHCOPYWARNING: AbstractContextManager = pd.option_context(
        'mode.chained_assignment', None
    )

    # Get-only properties

    @property
    def template(self) -> Dict[str, Tuple[str, ...]]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"template"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].template

    @property
    def mol_type(self) -> str:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"mol_type"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].mol_type

    @property
    def description(self) -> str:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"description"``]."""
        return self._WORKFLOW_TEMPLATES[self.name].description

    @property
    def import_columns(self) -> Dict[str, Tuple[str, str]]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"import_columns"``]."""  # noqa
        return self._WORKFLOW_TEMPLATES[self.name].import_columns

    @property
    def export_columns(self) -> Tuple[Tuple[str, str], ...]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"export_columns"``]."""  # noqa
        return self._WORKFLOW_TEMPLATES[self.name].export_columns

    # Getter and setter properties

    @property
    def read(self) -> bool:
        """Get or set :attr:`WorkFlow.read`.

        Setting accepts either a boolean or a container that may
        or may not contain :attr:`WorkFlow.mol_type` as value.

        """
        return self._read

    @read.setter
    def read(self, value: Union[bool, Container]) -> None:
        try:
            self._read = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._read = bool(value)

    @property
    def write(self) -> bool:
        """Get or set :attr:`WorkFlow.write`.

        Setting accepts either a boolean or a container that may
        or may not contain :attr:`WorkFlow.mol_type` as value.

        """
        return self._write

    @write.setter
    def write(self, value: Union[bool, Container]) -> None:
        try:
            self._write = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._write = bool(value)

    @property
    def overwrite(self) -> bool:
        """Get or set :attr:`WorkFlow.overwrite`.

        Setting accepts either a boolean or a container that may
        or may not contain :attr:`WorkFlow.mol_type` as value.

        """
        return self._overwrite

    @overwrite.setter
    def overwrite(self, value: Union[bool, Container]) -> None:
        try:
            self._overwrite = bool(self.db) and self.mol_type in value
        except TypeError:  # value is not a container
            self._overwrite = bool(value)

    @property
    def jobs(self) -> Tuple[Optional[Type[Job]], ...]:
        """Get or set :attr:`WorkFlow.jobs`.

        Setting accepts either a |plams.Job| type, ``None`` or
        an iterable containing one (or both) of the aforementioned objects.

        """
        return self._jobs

    @jobs.setter
    def jobs(self, value: Union[None, Type[Job], Iterable[None], Iterable[Type[Job]]]) -> None:
        if isinstance(value, Job):
            self._jobs = (value,)
        else:
            self._jobs = (None,) if value is None else tuple(value)

    @property
    def settings(self) -> Tuple[Optional[Settings], ...]:
        """Get or set :attr:`WorkFlow.settings`.

        Setting accepts either a |plams.Settings| instance, ``None`` or
        an iterable containing one (or both) of the aforementioned objects.

        """
        return self._settings

    @settings.setter
    def settings(self, value: Union[None, Settings, Iterable[None], Iterable[Settings]]) -> None:
        if isinstance(value, Settings):
            self._settings = (value,)
        else:
            self._settings = (None,) if value is None else tuple(value)

    # Methods and magic methods

    def __init__(self, name: str,
                 db: Optional['Database'] = None,
                 read: Union[bool, Container] = False,
                 write: Union[bool, Container] = False,
                 overwrite: Union[bool, Container] = False,
                 path: Optional[str] = None,
                 keep_files: bool = True,
                 read_template: bool = True,
                 jobs: OptionalJobType = None,
                 settings: OptionalSettings = None,
                 **kwargs: Any) -> None:
        """Initialize a :class:`WorkFlow` instance; see also :meth:`Workflow.from_template`."""
        super().__init__()

        if name not in self._WORKFLOW_TEMPLATES:
            err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                   f"Allowed values: {', '.join(repr(k) for k in self._WORKFLOW_TEMPLATES)}")
            raise ValueError(err)

        self.name: str = name
        self.db = db

        self.read: bool = read
        self.write: bool = write
        self.overwrite: bool = overwrite

        self.path: str = path if path is not None else os.getcwd()
        self.keep_files: bool = keep_files
        self.read_template: bool = read_template
        self.jobs: Iterable[Job] = jobs
        self.settings: Iterable[Settings] = settings

        for k, v in kwargs.items():
            if hasattr(self, k):
                raise AttributeError(f"An attribute by the name of '{k}' already exists")
            setattr(self, k, v)

    @AbstractDataClass.inherit_annotations()
    def _str_iterator(self):
        iterator = super()._str_iterator()
        return ((k.strip('_'), v) for k, v in iterator)

    def __call__(self, func: Callable, df: pd.DataFrame,
                 index: Union[slice, pd.Series] = slice(None),
                 columns: Optional[List[Hashable]] = None, **kwargs) -> None:
        r"""Initialize the workflow.

        Parameters
        ----------
        func : :data:`Callable<typing.Callable>`
            A callable object which will recieve **df**, all :class:`WorkFlow` instance
            attributes and ***kwargs** as arguments.
            The callable is expected to conduct some operation and export the results **dfs'**
            :attr:`WorkFlow.import_columns` columns.

        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        index : :class:`slice` or :class:`pandas.Series` [:class:`bool`]
            An object for slicing the rows of **df** (*i.e.* a :attr:`pandas.DataFrame.index`).

        columns : :class:`list` [:data:`Hashable<typing.Hashable>`], optional
            An object for slicing the columns of **df** (*i.e.* :attr:`pandas.DataFrame.columns`).
            The output of **func** will be fed into :code:`df[columns]`.
            If ``None``, use :attr:`WorkFlow.import_columns` instead.

        \**kwargs : :data:`Any<typing.Any>`
            Optional keyword arguments for **func**.

        See Also
        --------
        :meth:`Workflow.from_db`:
            Returns a value for the **idx_slice** parameter.

        """
        # Prepare slices
        slice1 = index, MOL
        slice2 = index, list(self.import_columns.keys()) if columns is None else columns

        # Run the workflow
        logger.info(f"Starting {self.description}")
        with PlamsInit(path=self.path, folder=self.name), self._SUPRESS_SETTINGWITHCOPYWARNING:
            self_vars = {k.strip('_'): v for k, v in vars(self).items()}
            value = func(df.loc[slice1], **self_vars, **kwargs)
            df.loc[slice2] = value
        logger.info(f"Finishing {self.description}\n")

    def from_db(self, df: pd.DataFrame, inplace: bool = True, get_mol: bool = True,
                columns: Optional[Dict[Hashable, Any]] = None) -> Union[slice, pd.Series]:
        """Ensure that all required keys are present in **df** and import from the database.

        Returns a :class:`pandas.index` with all to-be updated rows, as based on how many
        previous results were imported from :attr:`WorkFlow.db`.
        If no results were pulled from :attr:`WorkFlow.db` (or :attr:`WorkFlow.overwrite` is
        ``True``), then return :code:`slice(None)`.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        inplace : :class:`bool`
            If ``True``, perform an inplace update of the Cartesian coordinates of all molecules
            rather than importing new molecules.

        get_mol : :class:`bool`
            If ``False`` do *not* try to import molecules from the database.

        columns : :class:`dict` [:data:`Hashable<typing.Hashable>`, :class:`object`], optional
            An dictionary whose keys will be used for slicing the columns of **df**
            (*i.e.* :attr:`pandas.DataFrame.columns`).
            The dictionary valeus are used as fill values if a key belongs
            to a to-be created column.
            Recommended values are :data:`numpy.nan`, ``None``, ``-1`` and/or ``False``.
            If ``None``, use :attr:`WorkFlow.import_columns` instead.

        Returns
        -------
        :class:`pandas.Series` [:class:`bool`] or :class:`slice`
            A Series for slicing a part of **df** or  a :class:`slice` object for
            slicing the entirety of **df** (*i.e.* :code:`slice(0, None`).

        """
        # Add all necasary keys to **df**
        import_columns = self.import_columns if columns is None else columns
        for key, value in import_columns.items():
            if key not in df:
                df[key] = value

        if not self.read:  # Nothing to see here, move along
            return slice(None)

        # Import from the database
        with self._SUPRESS_SETTINGWITHCOPYWARNING:
            mol_list = self.db.from_csv(df, database=self.mol_type,
                                        get_mol=get_mol, inplace=inplace)
            if not inplace:  # mol_list is an actual sequence instead of None
                df[MOL] = mol_list

        # Return a new DataFrame slice based on previously calculated results
        if self.overwrite:
            return slice(None)
        else:
            keys = list(import_columns.keys())
            return self._isnull(df, keys).any(axis=1)

    def to_db(self, df: pd.DataFrame, status: Optional[str] = None,
              job_recipe: Optional[dict] = None,
              index: Union[slice, pd.Series] = slice(None),
              columns: Optional[Dict[Hashable, Any]] = None) -> None:
        """Export results to the database.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.

        status : :class:`str`, optional
            Whether or not **df** contains structures resulting from a geometry optimization.

        job_recipe : :class:`dict`
            A (nested) dictionary with the used job settings.

        index : :class:`slice` or :class:`pandas.Series` [:class:`bool`]
            An object for slicing the rows of **df** (*i.e.* :attr:`pandas.DataFrame.index`).

        columns : :class:`list` [:data:`Hashable<typing.Hashable>`], optional
            An object for slicing the columns of **df** (*i.e.* :attr:`pandas.DataFrame.columns`).
            If ``None``, use :attr:`WorkFlow.export_columns` instead.

        """
        # Dont export any settings columns if job_recipe is None
        # No job recipe == no settings to export anyway
        _export_columns = self.export_columns if columns is None else columns
        if job_recipe is None:
            export_columns = [i for i in _export_columns if i[0] != 'settings']
        else:
            export_columns = list(_export_columns)

        # Set the optimization status of the molecules to True
        if status == 'optimized':
            df.loc[index, OPT] = False

        # Write results to the database
        if self.write:
            with self._SUPRESS_SETTINGWITHCOPYWARNING:
                self.db.update_csv(
                    df.loc[index],
                    database=self.mol_type,
                    columns=export_columns,
                    overwrite=self.overwrite,
                    job_recipe=job_recipe,
                    status=status,
                )

        # Remove the PLAMS results directories
        if not self.keep_files:
            rmtree(Path(self.path) / self.name)

    @classmethod
    def from_template(cls, settings: Union[Settings, SettingsDataFrame], name: str) -> 'WorkFlow':
        """Construct a :class:`WorkFlow` instance from a |plams.Settings| object.

        Parameters
        ----------
        settings : |plams.Settings|
            A Settings instance with all CAT settings.
            Certain values are extracted from **settings** based on the supplied template
            (see **name**).

        name : :class:`str`
            The name of the settings template.

        See Also
        --------
        :attr:`WorkFlow._WORKFLOW_TEMPLATES`
            A dictionary with all available template names (*i.e.* its keys).

        """
        # Extract the settings object from the SettingsDataFrame
        if isinstance(settings, SettingsDataFrame):
            settings = settings.settings

        kwargs = {'name': name}

        # Raise a KeyError if a key cannot be found
        with Settings.supress_missing():
            try:  # Extract the correct template
                template: Dict[str, Tuple[str, ...]] = cls._WORKFLOW_TEMPLATES[name].template
            except KeyError as ex:
                err = (f"Invalid value for the 'name' parameter: {repr(name)}\n"
                       f"Allowed values: {', '.join(repr(k) for k in cls._WORKFLOW_TEMPLATES)}")
                raise ValueError(err).with_traceback(ex.__traceback__)

            # Create a dictionary with keyword arguments
            for k, v in template.items():
                kwargs[k] = settings.get_nested(v)

        # Post process all jobs and job settings
        kwargs['jobs'] = pop_and_concatenate(kwargs, 'job')
        kwargs['settings'] = pop_and_concatenate(kwargs, 's')
        return cls.from_dict(kwargs)

    def get_recipe(self) -> Settings:
        """Create a recipe for :meth:`WorkFlow.to_db`."""
        settings_names = [i[1:] for i in self.export_columns if i[0] == 'settings']
        uff_fallback = {
            'key': f'RDKit_{rdkit.__version__}', 'value': f'{UFF.__module__}.{UFF.__name__}'
        }

        ret = Settings()
        for name, job, settings in zip(settings_names, self.jobs, self.settings):
            # job is None, *i.e.* it's an RDKit UFF optimziation
            if job is None:
                ret[name].update(uff_fallback)
                continue

            settings = Settings(settings)
            if self.read_template:  # Update the settings using a QMFlows template
                template = qmflows.geometry['specific'][self.type_to_string(job)].copy()
                settings.soft_update(template)
            ret[name].key = job
            ret[name].value = settings
        return ret

    @staticmethod
    def _isnull(df: pd.DataFrame, columns: List[Hashable]) -> pd.DataFrame:
        """A more expansive version of the :func:`pandas.isnull` function.

        :class:`int` series now also return ``True`` if smaller than ``0`` and :class:`bool`
        series are simply inverted.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame.

        columns : :class:`list`
            A list of column keys from **df**.

        """
        dtype_dict: Dict[np.dtype, Callable] = {
            np.dtype(bool): operator.invert,
            np.dtype(int): _lt_0,
            np.dtype(float): pd.isnull,
            np.dtype(object): pd.isnull
        }

        ret = pd.DataFrame(index=df.index)
        for key, series in df[columns].items():
            func = dtype_dict.get(series.dtype, pd.isnull)
            ret[key] = func(series)
        return ret

    @staticmethod
    def type_to_string(job: Union[Job, Type[Job]]) -> Optional[None]:
        """Turn a :class:`type` instance into a :class:`str`.

        Accepts one of the following |plams.Job| subclasses:
            * ``ADFJob``
            * ``AMSJob``
            * ``DiracJob``
            * ``Cp2kJob``
            * ``GamessJob``
            * ``ORCAJob``

        Parameters
        ----------
        job : :class:`type` [|plams.Job|] or |plams.Job|
            A PLAMS Job type or instance.

        Returns
        -------
        :class:`str`, optional
            Returns either ``None`` or an item pulled from :data:`._job_dict`.

        """
        if not isinstance(job, type):
            job = type(job)  # Convert a class instance into a class

        try:
            return _job_dict[job]
        except KeyError as ex:
            logger.error(f"No default settings available for type: '{job.__class__.__name__}'")
            logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)
            return None

    @staticmethod
    def pop_job_settings(mol_list: Iterable[Molecule], key: str = 'job_path') -> List[List[str]]:
        """Take a list of molecules and pop and return all references to **key**.

        Parameters
        ----------
        mol_list : :data:`Iterable<typing.Iterable>` [|plams.Molecule|]
            An iterable consisting of PLAMS molecules.
            For this method to be effective they should contain a property by the name of **key**:
            a list of strings represnting paths to .in files.

        key : :data:`Hashable<typing.Hashable>`
            The to-be popped key in each molecule in **mol_list**.

        Returns
        -------
        :class:`list` [:class:`list` [:class:`str`]]
            A nested list of strings popped from **mol_list**.

        """
        ret = []
        for mol in mol_list:
            try:
                ret.append(mol.properties.pop(key))
            except KeyError:
                ret.append([])
            mol.properties[key] = []
        return ret


class PlamsInit(AbstractContextManager):
    """A context manager for calling :func:`.restart_init` and |plams.finish|."""

    def __init__(self, path: str, folder: str, hashing: str = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        """Enter the context manager; call :func:`.restart_init`."""
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; call |plams.finish|."""
        finish()