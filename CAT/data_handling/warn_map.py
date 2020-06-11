"""A module for managing exceptions and warnings in CAT."""

import warnings
from types import MappingProxyType
from typing import Mapping, Callable, NoReturn, Union, Type

from scm.plams import MoleculeError

__all__ = ['WARN_MAP']


class MoleculeWarning(Warning, MoleculeError):
    """A :exc:`Warning` subclass for :class:`Molecule<scm.plams.mol.molecule.Molecule>` related errors."""  # noqa


class ValueWarning(Warning, ValueError):
    """A :exc:`Warning` subclass for :exc:`ValueError` related errors."""  # noqa


#: Map an :exc:`Exception` type to a :exc:`Warning` type.
CATEGORY_MAP: Mapping[Type[Exception], Type[Warning]] = MappingProxyType({
    MoleculeError: MoleculeWarning,
    ValueError: ValueWarning
})


def _warn(exc: Exception) -> None:
    """Perform a warning using **exc**.

    When possible, the warning category will be derived from the passed Exception type
    (see :data:`CATEGORY_MAP`).
    Will default to :exc:`RuntimeWarning` otherwise.

    """  # noqa
    warnings.warn(str(exc), stacklevel=2,
                  category=CATEGORY_MAP.get(type(exc), RuntimeWarning))


def _raise(exc: Exception) -> NoReturn:
    """Raise **exc**."""
    raise exc


def _ignore(exc: Exception) -> None:
    """Do nothing."""
    return None


#: Map a string to callable for either raising an :exc:`Exception`,
#: displaying a :exc:`Warning` or doing nothing.
WARN_MAP: Mapping[str, Callable[[Exception], Union[NoReturn, None]]] = MappingProxyType({
    'raise': _raise,
    'warn': _warn,
    'ignore': _ignore
})
