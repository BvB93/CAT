"""A module which adds the :class:`.FrozenSettings` class, an immutable counterpart to plams.Settings_.

.. _plams.Settings: https://www.scm.com/doc/plams/components/settings.html

Index
-----
.. currentmodule:: CAT.frozen_settings
.. autosummary::
    FrozenSettings

API
---
.. autoclass:: FrozenSettings
    :members:
    :private-members:
    :special-members:

"""  # noqa: E501

import copy
import textwrap
from typing import (NoReturn, Hashable)

from scm.plams import Settings

__all__ = ['FrozenSettings']


class FrozenSettings(Settings):
    """An inmutable subclass of plams.Settings_.

    .. _plams.Settings: https://www.scm.com/doc/plams/components/settings.html
    """

    def __init__(self, *args, **kwargs) -> None:
        """Initialize the construction of a :class:`FrozenSettings` instance."""
        dict.__init__(self, *args, **kwargs)
        cls = type(self)

        # Fill the FrozenSettings instance by means of the dict.__setitem__ method
        for key, value in self.items():
            if isinstance(value, dict):
                value_new = cls(value)
            elif isinstance(value, list):
                value_new = [cls(i) if isinstance(i, dict) else i for i in value]
            else:
                continue
            dict.__setitem__(self, key, value_new)

        # An attribute for caching the hash of this instance
        dict.__setattr__(self, '_hash', None)

    def __str__(self) -> str:
        """Return a string representation of this instance."""
        if not self:
            return f'{self.__class__.__name__}()'

        indent = 4 * ' '
        ret = super().__str__()[:-1]
        return f'{self.__class__.__name__}(\n{textwrap.indent(ret, indent)}\n)'

    __repr__ = __str__

    def __missing__(self, key: Hashable) -> 'FrozenSettings':
        """Return :data:`_frozen_settings`."""
        return _frozen_settings

    def __bool__(self) -> bool:
        """Implement :class:`bool(self)<bool>`."""
        return bool(len(self))

    @classmethod
    def _raise_exc(cls, *args, **kwargs) -> NoReturn:
        """Method unsupported; raise a :exc:`TypeError`."""
        raise TypeError(f"'{cls.__name__}' object does not support item assignment or deletion")

    __delitem__ = __setitem__ = set_nested = _raise_exc

    def flatten(self, flatten_list: bool = True) -> 'FrozenSettings':
        """Return a flattened copy of this instance.

        New keys are constructed by concatenating the (nested) keys of this instance into tuples.
        Opposite of the :meth:`.Settings.unflatten` method.
        If *flatten_list* is ``True``, all nested lists will be flattened as well. Dictionary keys are replaced with list indices in such case.

        """  # noqa: E501
        ret = super().flatten(flatten_list)
        return type(self)(ret)

    def unflatten(self, unflatten_list: bool = True) -> 'FrozenSettings':
        """Return a nested copy of this instance.

        New keys are constructed by expanding the keys of this instance (*e.g.* tuples) into new nested |Settings| instances.
        If *unflatten_list* is ``True``, integers will be interpretted as list indices and are used for creating nested lists.
        Opposite of the :meth:`.flatten` method.
        """  # noqa: E501
        ret = super().unflatten(unflatten_list)
        return type(self)(ret)

    def __hash__(self) -> int:
        """Return the hash of this instance; pull it from :attr:`._hash` if possible."""
        if self._hash is not None:
            return dict.__getattribute__(self, '_hash')

        flat = super().flatten()
        ret = 0
        for k, v in flat.items():
            ret ^= hash(k + (v, ))
        dict.__setattr__(self, '_hash', ret)
        return ret

    def copy(self, deep: bool = False) -> 'FrozenSettings':
        """Create a copy of this instance.

        The returned instance is a recursive copy, the **deep** parameter only affecting the
        function used for copying non-:class:`dict` values:
        :func:`copy.copy` if ``deep=False`` and  :func:`copy.deepcopy` if ``deep=True``.

        """
        copy_func = copy.deepcopy if deep else copy.copy
        ret = type(self)()
        for key, value in self.items():
            value_cp = copy_func(value)
            Settings.__setitem__(ret, key, value_cp)
        return ret

    def __copy__(self) -> 'FrozenSettings':
        """Create a shallow copy of this instance by calling :meth:`FrozenSettings.copy`."""
        return self.copy(deep=False)

    def __deepcopy__(self, memo=None) -> 'FrozenSettings':
        """Create a deep copy of this instance by calling :meth:`FrozenSettings.copy`."""
        return self.copy(deep=True)


_frozen_settings = FrozenSettings()
