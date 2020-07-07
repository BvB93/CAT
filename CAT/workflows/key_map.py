"""A module for storing all column keys used throughout CAT.

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

#: A :class:`~collections.abc.Mapping` containing all column keys used throughout CAT.
KEY_MAP: Mapping[str, Tuple[str, str]] = MappingProxyType({
    'OPT': ('opt', ''),
    'MOL': ('mol', ''),
    'FORMULA': ('formula', ''),
    'HDF5_INDEX': ('hdf5 index', ''),
    'JOB_SETTINGS_QD_OPT': ('settings', 'qd_opt'),
    'JOB_SETTINGS_CRS': ('settings', 'CRS'),
    'JOB_SETTINGS_BDE': ('settings', 'BDE'),
    'JOB_SETTINGS_ASA': ('settings', 'ASA'),
    'JOB_SETTINGS_CDFT': ('settings', 'CDFT'),
    'CDFT_MU': ('cdft', 'Electronic chemical potential (mu)'),
    'CDFT_CHI': ('cdft', 'Electronegativity (chi=-mu)'),
    'CDFT_ETA': ('cdft', 'Hardness (eta)'),
    'CDFT_S': ('cdft', 'Softness (S)'),
    'CDFT_GAMMA': ('cdft', 'Hyperhardness (gamma)'),
    'CDFT_OMEGA': ('cdft', 'Electrophilicity index (w=omega)'),
    'CDFT_NUCLEOFUGE': ('cdft', 'Dissocation energy (nucleofuge)'),
    'CDFT_ELECTROFUGE': ('cdft', 'Dissociation energy (electrofuge)'),
    'CDFT_W_MINUS': ('cdft', 'Electrodonating power (w-)'),
    'CDFT_W_PLUS': ('cdft', 'Electroaccepting power(w+)'),
    'CDFT_ELECTROPHILICITY': ('cdft', 'Net Electrophilicity'),
    'CDFT_DELTAF_MINUS': ('cdft', 'Global Dual Descriptor Deltaf-'),
    'CDFT_DELTAF_PLUS': ('cdft', 'Global Dual Descriptor Deltaf+'),
    'CDFT_MU_MINUS': ('cdft', 'Electronic chemical potential (mu-)'),
    'CDFT_MU_PLUS': ('cdft', 'Electronic chemical potential (mu+)'),
    'ASA_INT': ('ASA', 'E_int'),
    'ASA_STRAIN': ('ASA', 'E_strain'),
    'ASA_E': ('ASA', 'E'),
    'V_BULK': ('V_bulk', ''),
    'LIGAND_COUNT': ('ligand count', '')
})

globals().update(KEY_MAP)

__all__ = list(KEY_MAP.keys())
