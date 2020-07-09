"""A module for holding the :class:`WorkFlow` class.

Index
-----
.. currentmodule:: CAT.workflows.workflow
.. autosummary::
    WorkFlow

API
---
.. autoclass:: WorkFlow
    :members:

"""

import os
from shutil import rmtree
from pathlib import Path
from collections import abc
from typing import (
    Optional, Union, Dict, MutableMapping, TypeVar, Iterable, Container, Tuple, Callable,
    Any, List, Type, Mapping, TYPE_CHECKING, cast, ClassVar, ContextManager
)

import numpy as np
import pandas as pd

from scm.plams import finish, Settings, Molecule
from scm.plams.core.basejob import Job
from assertionlib import AbstractDataClass, NDRepr
from nanoutils import Literal

from .key_map import MOL, OPT, HDF5_INDEX
from .workflow_dicts import WORKFLOW_TEMPLATE, _TemplateMapping
from ..utils import restart_init
from ..logger import logger
from ..settings_dataframe import SettingsDataFrame

if TYPE_CHECKING:
    from dataCAT import Database
    from numpy.typing import ArrayLike
else:
    Database = 'dataCAT.Database'
    ArrayLike = 'numpy.typing.ArrayLike'

__all__ = ['WorkFlow']

NDRepr.repr_SettingsDataFrame = NDRepr.repr_DataFrame  # type: ignore
aNDRepr = NDRepr()

T = TypeVar('T')


def _return_true(value: object) -> bool:
    """Return :data:`True`."""
    return True


def _lt_0(value) -> bool:
    """Return if **value** is smaller than ``0``."""
    return value < 0


def pop_and_concatenate(mapping: MutableMapping[str, T], base_key: object,
                        filter_func: Callable[[Any], bool] = _return_true) -> Tuple[T, ...]:
    """Take a key and :meth:`pop<dict.pop>` all values from **mapping**.

    The popping will continue as long as :code:`base_key + str(i)` is available in the mapping,
    where ``i`` is an :class:`int` larger than 1.
    The value if ``i`` will start from 1 and increase by `+ 1` every iteration.

    Examples
    --------
    .. code:: python

        >>> from CAT.workflows.workflow import pop_and_concatenate

        >>> mapping: dict = {
        ...     'job1': 1,
        ...     'job2': 2,
        ...     'job3': 3,
        ...     'final_key': True
        ... }

        >>> base_key = 'job'
        >>> value_tuple: tuple = pop_and_concatenate(mapping, base_key)
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
    ret: List[T] = []
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


OptionalJobType = Union[None, Type[Job], Iterable[Optional[Type[Job]]]]
OptionalSettings = Union[None, Settings, Iterable[Optional[Settings]]]


class WorkFlow(AbstractDataClass):
    """A workflow manager.

    Examples
    --------
    Typical usage example:

    .. code:: python

        >>> import pandas as pd

        # Prepare workflow parameters
        >>> df = pd.DataFrame(...)  # doctest: +SKIP
        >>> settings = Settings()

        >>> def fancy_df_func(df, **kwargs):
        ...     pass

        # Create the workflow
        >>> workflow = WorkFlow.from_template(settings, name='asa')  # doctest: +SKIP
        >>> print(workflow)  # doctest: +SKIP
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

        # Run the workflow
        >>> idx = workflow.from_db(df)  # doctest: +SKIP
        >>> workflow(fancy_df_func, df, index=idx)  # doctest: +SKIP

        # Export all workflow results
        >>> job_recipe = workflow.get_recipe()  # doctest: +SKIP
        >>> workflow.to_db(df, job_recip=job_recipe)  # doctest: +SKIP

    """

    #: Map a name to a workflow template.
    _WORKFLOW_TEMPLATES: ClassVar[Mapping[str, _TemplateMapping]] = WORKFLOW_TEMPLATE

    #: A context manager for supressing Pandas :exc:`SettingwithCopyWarning`.
    _SUPRESS_SETTINGWITHCOPYWARNING: ClassVar[ContextManager[None]] = pd.option_context(
        'mode.chained_assignment', None
    )

    # Get-only properties

    @property
    def template(self) -> Mapping[str, Tuple[str, ...]]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"template"``]."""
        return self._WORKFLOW_TEMPLATES[self.name]['template']

    @property
    def mol_type(self) -> Literal['core', 'ligand', 'qd']:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"mol_type"``]."""
        return self._WORKFLOW_TEMPLATES[self.name]['mol_type']

    @property
    def description(self) -> str:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"description"``]."""
        return self._WORKFLOW_TEMPLATES[self.name]['description']

    @property
    def import_columns(self) -> Mapping[Tuple[str, str], np.generic]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"import_columns"``]."""  # noqa: E501
        return self._WORKFLOW_TEMPLATES[self.name]['import_columns']

    @property
    def export_columns(self) -> Tuple[Tuple[str, str], ...]:
        """Get :attr:`WorkFlow._WORKFLOW_TEMPLATES` [:attr:`WorkFlow.name`] [``"export_columns"``]."""  # noqa: E501
        return self._WORKFLOW_TEMPLATES[self.name]['export_columns']

    # Getter and setter properties

    @property
    def read(self) -> bool:
        """Get or set :attr:`WorkFlow.read`.

        Setting accepts either a boolean or a container that may
        or may not contain :attr:`WorkFlow.mol_type` as value.

        """
        return self._read

    @read.setter
    def read(self, value: Union[bool, Container[str]]) -> None:
        try:
            self._read: bool = bool(self.db) and self.mol_type in value  # type: ignore
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
    def write(self, value: Union[bool, Container[str]]) -> None:
        try:
            self._write: bool = bool(self.db) and self.mol_type in value
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
    def overwrite(self, value: Union[bool, Container[str]]) -> None:
        try:
            self._overwrite: bool = bool(self.db) and self.mol_type in value
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
    def jobs(self, value: Union[None, Type[Job], Iterable[Optional[Type[Job]]]]) -> None:
        if isinstance(value, type):
            self._jobs: Tuple[Optional[Type[Job]], ...] = (value,)
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
    def settings(self, value: Union[None, Settings, Iterable[Optional[Settings]]]) -> None:
        if isinstance(value, Settings):
            self._settings: Tuple[Optional[Settings], ...] = (value,)
        else:
            self._settings = (None,) if value is None else tuple(value)

    # Methods and magic methods

    def __init__(self, name: str,
                 db: Optional[Database] = None,
                 read: Union[bool, Container[str]] = False,
                 write: Union[bool, Container[str]] = False,
                 overwrite: Union[bool, Container[str]] = False,
                 path: Union[None, str, 'os.PathLike[str]'] = None,
                 keep_files: bool = True,
                 jobs: OptionalJobType = None,
                 settings: OptionalSettings = None,
                 **kwargs: Any) -> None:
        """Initialize a :class:`WorkFlow` instance; see also :meth:`Workflow.from_template`."""
        super().__init__()

        if name not in self._WORKFLOW_TEMPLATES:
            raise ValueError(f"Invalid value for the 'name' parameter: {name!r}\n"
                             f"Allowed values: {list(self._WORKFLOW_TEMPLATES.keys())!r}")

        self.name = name
        self.db = db

        self.read = cast(bool, read)
        self.write = cast(bool, write)
        self.overwrite = cast(bool, overwrite)

        self.path: Union[str, 'os.PathLike[str]'] = path if path is not None else os.getcwd()
        self.keep_files = keep_files
        self.jobs = cast(Tuple[Optional[Type[Job]], ...], jobs)
        self.settings = cast(Tuple[Optional[Settings], ...], settings)

        for k, v in kwargs.items():
            if hasattr(self, k):
                raise AttributeError(f"An attribute by the name of {k!r} already exists")
            setattr(self, k, v)

    @AbstractDataClass.inherit_annotations()
    def _str_iterator(self):
        iterator = super()._str_iterator()
        return ((k.strip('_'), v) for k, v in iterator)

    # TODO: Ensure that func allways returns an array-like object (no iterators)
    def __call__(self, func: Callable, df: pd.DataFrame,
                 index: Union[slice, pd.Series] = slice(None),
                 columns: Optional[ArrayLike] = None,
                 no_loc: bool = False, **kwargs) -> None:
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

        no_loc : :class:`bool`
            If ``True``, substitute :meth:`pandas.DataFrame.loc` for
            :meth:`pandas.DataFrame.__setitem__`.
            Usefull to prevent Panda's rather aggressive typecasting when passing object arrays.
            Note that **no_loc** and **index** cannot be simultaneously specified.

        \**kwargs : :data:`Any<typing.Any>`
            Optional keyword arguments for **func**.

        See Also
        --------
        :meth:`Workflow.from_db`:
            Returns a value for the **idx_slice** parameter.

        """
        # Prepare slices
        if no_loc and index != slice(None):
            raise ValueError
        slice1 = index, MOL
        slice2 = index, list(self.import_columns.keys()) if columns is None else columns

        # Run the workflow
        logger.info(f"Starting {self.description}")
        with PlamsInit(path=self.path, folder=self.name), self._SUPRESS_SETTINGWITHCOPYWARNING:
            self_vars = {k.strip('_'): v for k, v in vars(self).items()}
            value = func(df.loc[slice1], columns=columns, **self_vars, **kwargs)

            if not isinstance(value, abc.Iterator) and not np.any(value):
                return
            elif no_loc:
                for k, v in zip(slice2[1], value):
                    df[k] = v
            else:
                try:
                    df.loc[slice2] = value
                except ValueError as ex:
                    logger.debug(f"df = {aNDRepr.repr(df)}")
                    logger.debug(f"index = {aNDRepr.repr(slice2[0])}")
                    logger.debug(f"columns = {aNDRepr.repr(slice2[1])}")
                    logger.debug(f"value = {aNDRepr.repr(value)}")
                    raise ex
        logger.info(f"Finishing {self.description}\n")

    def from_db(self, df: pd.DataFrame, *columns: str, read_mol: bool = False) -> pd.DataFrame:
        """Ensure that all required keys are present in **df** and import from the database.

        Returns a :class:`pandas.index` with all to-be updated rows, as based on how many
        previous results were imported from :attr:`WorkFlow.db`.
        If no results were pulled from :attr:`WorkFlow.db` (or :attr:`WorkFlow.overwrite` is
        ``True``), then return :code:`slice(None)`.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A DataFrame with molecules and results.
        get_mol : :class:`bool`
            If ``False`` do *not* try to import molecules from the database.

        Returns
        -------
        :class:`pandas.Series[bool]<pandas.Series>` or :class:`slice`
            A Series for slicing a part of **df** or  a :class:`slice` object for
            slicing the entirety of **df** (*i.e.* :code:`slice(0, None`).

        """
        if self.db is None:
            df_bool = df.copy()
            df_bool[:] = True
            return df_bool

        if not self.read:
            _df = df
            df = df.index

        df, df_bool = self.db.to_df(df, self.mol_type, *columns, read_mol=read_mol)

        if self.overwrite:
            df_bool[:] = True
        if not self.read and read_mol:
            _df[HDF5_INDEX] = df[HDF5_INDEX]
        return df_bool

    def to_db(self, df: pd.DataFrame, df_bool: pd.DataFrame,
              columns: Optional[ArrayLike] = None,
              status: Optional[str] = None) -> None:
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
        if self.db is None:
            return

        # Set the optimization status of the molecules to True
        if status == 'optimized':
            df[OPT] = True

        # The group name
        name = self.mol_type
        if status == 'no_opt':
            status = None
            name += '_no_opt'  # type: ignore

        # Write results to the database
        if self.write:
            with self._SUPRESS_SETTINGWITHCOPYWARNING:
                self.db.from_df(
                    df, df_bool, name,
                    columns=columns,
                    overwrite=self.overwrite,
                    status=status
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

        kwargs: Dict[str, Any] = {'name': name}

        # Raise a KeyError if a key cannot be found
        with Settings.supress_missing():
            try:  # Extract the correct template
                template = cls._WORKFLOW_TEMPLATES[name]['template']
            except KeyError as ex:
                err = (f"Invalid value for the 'name' parameter: {name!r}\n"
                       f"Allowed values: {', '.join(repr(k) for k in cls._WORKFLOW_TEMPLATES)}")
                raise ValueError(err) from ex

            # Create a dictionary with keyword arguments
            for k, v in template.items():
                kwargs[k] = settings.get_nested(v)

        # Post process all jobs and job settings
        kwargs['jobs'] = pop_and_concatenate(kwargs, 'job')
        kwargs['settings'] = pop_and_concatenate(kwargs, 's')
        return cls.from_dict(kwargs)

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
            ret.append(mol.properties.pop(key, []))
            mol.properties[key] = []
        return ret


class PlamsInit(ContextManager[None]):
    """A context manager for calling :func:`.restart_init` and |plams.finish|."""

    def __init__(self, path: Union[str, 'os.PathLike[str]'],
                 folder: Union[str, 'os.PathLike[str]'],
                 hashing: str = 'input'):
        self.path = path
        self.folder = folder
        self.hashing = hashing

    def __enter__(self) -> None:
        """Enter the context manager; call :func:`.restart_init`."""
        restart_init(self.path, self.folder, self.hashing)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; call |plams.finish|."""
        finish()
