"""A module for holding the :class:`.SettingsDataFrame` and :class:`.SettingsSeries` classes.

Index
-----
.. currentmodule:: CAT.settings_dataframe
.. autosummary::
    SettingsSeries
    SettingsDataFrame

API
---
.. autoclass:: CAT.settings_dataframe.SettingsSeries
    :members:
    :private-members:
    :special-members:

.. autoclass:: CAT.settings_dataframe.SettingsDataFrame
    :members:
    :private-members:
    :special-members:

"""

from typing import Optional, Callable

import pandas as pd

from .frozen_settings import FrozenSettings

__all__ = ['SettingsSeries', 'SettingsDataFrame']


class SettingsSeries(pd.Series):
    """A subclass of the Pandas Series with an additional :attr:`.settings` attribute.

    Parameters
    ----------
    settings : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.SettingsSeries.settings`.

    Attributes
    ----------
    settings : |CAT.FrozenSettings|_
        An immutable :class:`.FrozenSettings` instance with additional user-defined (meta-)data.

    """

    _metadata = ['settings']

    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False,
                 settings: Optional[dict] = None) -> None:
        """Initialize the :class:`.SettingsSeries` construction."""
        self.settings = self._sanitize_settings(settings)
        super().__init__(data, index, dtype, name, copy, fastpath)

    @property
    def _constructor(self) -> Callable[..., 'SettingsSeries']:
        """Construct a :class:`.SettingsSeries` instance."""
        def _series(*args, **kwargs) -> SettingsSeries:
            return SettingsSeries(*args, **kwargs).__finalize__(self)
        return _series

    @property
    def _constructor_expanddim(self) -> Callable[..., 'SettingsDataFrame']:
        """Construct a :class:`.SettingsDataFrame` instance."""
        def _df(*args, **kwargs) -> SettingsDataFrame:
            return SettingsDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @staticmethod
    def _sanitize_settings(settings: Optional[dict]) -> FrozenSettings:
        """Sanitize the **settings** parameter for :attr:`SettingsSeries.settings`."""
        if settings is None:
            return FrozenSettings()
        elif isinstance(settings, dict):
            return FrozenSettings(settings)
        else:
            err = "The settings argument expects an instance of 'dict'; observed type '{}'"
            raise TypeError(err.format(settings.__class__.__name__))


class SettingsDataFrame(pd.DataFrame):
    """A subclass of the Pandas DataFrame with an additional :attr:`.settings` attribute.

    Parameters
    ----------
    settings : dict
        Optional: A dictionary with additional user-defined (meta-)data.
        See :attr:`.SettingsDataFrame.settings`.

    Attributes
    ----------
    settings : |CAT.FrozenSettings|_
        An immutable :class:`.FrozenSettings` instance with additional user-defined (meta-)data.

    """

    _metadata = ['settings']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False,
                 settings: Optional[dict] = None) -> None:
        """Initialize the :class:`.SettingsDataFrame` construction."""
        self.settings = self._sanitize_settings(settings)
        super().__init__(data, index, columns, dtype, copy)

    @property
    def _constructor(self) -> Callable[..., 'SettingsDataFrame']:
        """Construct a :class:`.SettingsDataFrame` instance."""
        def _df(*args, **kwargs) -> SettingsDataFrame:
            return SettingsDataFrame(*args, **kwargs).__finalize__(self)
        return _df

    @property
    def _constructor_sliced(self) -> Callable[..., 'SettingsSeries']:
        """Construct a :class:`.SettingsSeries` instance."""
        def _series(*args, **kwargs) -> SettingsSeries:
            return SettingsSeries(*args, **kwargs).__finalize__(self)
        return _series

    @staticmethod
    def _sanitize_settings(settings: Optional[dict]) -> FrozenSettings:
        """Sanitize the **settings** parameter for :attr:`SettingsDataFrame.settings`."""
        if settings is None:
            return FrozenSettings()
        elif isinstance(settings, dict):
            return FrozenSettings(settings)
        else:
            err = "The settings argument expects an instance of 'dict'; observed type '{}'"
            raise TypeError(err.format(settings.__class__.__name__))
