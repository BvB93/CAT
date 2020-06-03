"""
CAT.workflows.key_map
=====================

A module for storing all column keys used throughout CAT.

Index
-----
.. currentmodule:: CAT.workflows.key_map
.. autosummary::
    KEY_MAP

API
---
.. autodata:: KEY_MAP
    :annotation: : Mapping[str, Tuple[str, str]]

"""

from types import MappingProxyType
from typing import Mapping, Tuple

#: A :class:`Mapping<collections.abc.Mapping>` containing all column keys used throughout CAT.
KEY_MAP: Mapping[str, Tuple[str, str]] = MappingProxyType({
    'OPT': ('opt', ''),
    'MOL': ('mol', ''),
    'FORMULA': ('formula', ''),
    'HDF5_INDEX': ('hdf5 index', ''),
    'JOB_SETTINGS_QD_OPT': ('job_settings_qd_opt', ''),
    'JOB_SETTINGS_CRS': ('job_settings_crs', ''),
    'JOB_SETTINGS_BDE': ('job_settings_BDE', ''),
    'JOB_SETTINGS_ASA': ('job_settings_ASA', ''),
    'ASA_INT': ('ASA', 'E_int'),
    'ASA_STRAIN': ('ASA', 'E_strain'),
    'ASA_E': ('ASA', 'E'),
    'SETTINGS1': ('settings', '1'),
    'SETTINGS2': ('settings', '2'),
    'SETTINGS_SOLV1': ('settings', 'solv 1'),
    'SETTINGS_SOLV2': ('settings', 'solv 2'),
    'SETTINGS_ASA': ('settings', 'ASA 1'),
    'SETTINGS_BDE1': ('settings', 'BDE 1'),
    'SETTINGS_BDE2': ('settings', 'BDE 2'),
    'V_BULK': ('V_bulk', '')
})

globals().update(KEY_MAP)

__all__ = list(KEY_MAP.keys())
