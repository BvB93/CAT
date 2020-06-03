"""
CAT.workflows.workflow_dicts
============================

A module for loading :class:`WorkFlow` templates.

"""

import os
import sys
from types import MappingProxyType
from typing import Mapping, MutableMapping, Tuple, Dict, Any, TYPE_CHECKING

import yaml
import numpy as np

from .key_map import (
    OPT,
    HDF5_INDEX,
    JOB_SETTINGS_QD_OPT,
    JOB_SETTINGS_CRS,
    JOB_SETTINGS_BDE,
    JOB_SETTINGS_ASA,
    JOB_SETTINGS_CDFT,
    ASA_INT,
    ASA_STRAIN,
    ASA_E,
    SETTINGS1,
    SETTINGS2,
    SETTINGS_SOLV1,
    SETTINGS_SOLV2,
    SETTINGS_ASA,
    SETTINGS_BDE1,
    SETTINGS_BDE2,
    SETTINGS_CDFT,
    V_BULK,
    CDFT_MU,
    CDFT_CHI,
    CDFT_ETA,
    CDFT_S,
    CDFT_GAMMA,
    CDFT_OMEGA,
    CDFT_NUCLEOFUGE,
    CDFT_ELECTROFUGE,
    CDFT_W_MINUS,
    CDFT_W_PLUS,
    CDFT_ELECTROPHILICITY,
    CDFT_DELTAF_MINUS,
    CDFT_DELTAF_PLUS,
    CDFT_MU_MINUS,
    CDFT_MU_PLUS
)

if TYPE_CHECKING:
    from nanoutils import TypedDict

    class _TemplateMapping(TypedDict):
        description: str
        mol_type: str
        template: Mapping[str, Tuple[str, ...]]
        import_columns: Mapping[Tuple[str, str], float]
        export_columns: Tuple[Tuple[str, str], ...]

else:
    from typing import Mapping as _TemplateMapping

__all__ = ['WORKFLOW_TEMPLATE']


def _load_templates() -> Dict[str, Any]:
    """Load the templates from ``CAT/data/workflow_dicts/workflow_yaml.yaml``."""
    path = os.path.join(os.path.dirname(__file__), 'workflow_yaml.yaml')
    with open(path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


def _recursive_mapping_proxy(dct: MutableMapping) -> MappingProxyType:
    """Recursivelly convert **dct**, and all nested mutable mappings, into :class:`MappingProxyType<types.MappingProxyType>` instances."""  # noqa
    for k, v in dct.items():
        if isinstance(v, dict):
            dct[k] = _recursive_mapping_proxy(v)
    return MappingProxyType(dct)


def finalize_templates():
    """Update the templates using :data:`._TEMPLATE_UPDATE`."""
    base_dict = {
        'asa': {'import_columns': {ASA_INT: np.nan, ASA_STRAIN: np.nan, ASA_E: np.nan},
                'export_columns': (JOB_SETTINGS_ASA, SETTINGS_ASA, ASA_INT, ASA_STRAIN, ASA_E)},
        'ligand_opt': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                       'export_columns': (HDF5_INDEX, OPT, SETTINGS1, SETTINGS2)},
        'qd_attach': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                      'export_columns': (HDF5_INDEX,)},
        'qd_opt': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                   'export_columns': (HDF5_INDEX, OPT, JOB_SETTINGS_QD_OPT, SETTINGS1, SETTINGS2)},
        'crs': {'import_columns': {},
                'export_columns': (JOB_SETTINGS_CRS, SETTINGS_SOLV1, SETTINGS_SOLV2)},
        'bde': {'import_columns': {},
                'export_columns': (JOB_SETTINGS_BDE, SETTINGS_BDE1, SETTINGS_BDE2)},
        'forcefield': {'import_columns': {},
                       'export_columns': ()},
        'bulkiness': {'import_columns': {V_BULK: np.nan},
                      'export_columns': (V_BULK,)},
        'multi_ligand': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                         'export_columns': (HDF5_INDEX, OPT)},
        'cdft': {'import_columns': {CDFT_MU: np.nan, CDFT_CHI: np.nan, CDFT_ETA: np.nan,
                                    CDFT_S: np.nan, CDFT_GAMMA: np.nan, CDFT_OMEGA: np.nan,
                                    CDFT_NUCLEOFUGE: np.nan, CDFT_ELECTROFUGE: np.nan,
                                    CDFT_W_MINUS: np.nan, CDFT_W_PLUS: np.nan,
                                    CDFT_ELECTROPHILICITY: np.nan, CDFT_DELTAF_MINUS: np.nan,
                                    CDFT_DELTAF_PLUS: np.nan, CDFT_MU_MINUS: np.nan,
                                    CDFT_MU_PLUS: np.nan},
                 'export_columns': (JOB_SETTINGS_CDFT, SETTINGS_CDFT, CDFT_MU, CDFT_CHI, CDFT_ETA, CDFT_S,
                                    CDFT_GAMMA, CDFT_OMEGA, CDFT_NUCLEOFUGE, CDFT_ELECTROFUGE,
                                    CDFT_W_MINUS, CDFT_W_PLUS, CDFT_ELECTROPHILICITY,
                                    CDFT_DELTAF_MINUS, CDFT_DELTAF_PLUS, CDFT_MU_MINUS,
                                    CDFT_MU_PLUS)},
    }

    templates = _load_templates()
    for k, v1 in templates.items():
        v2 = base_dict[k]
        v1.update(v2)

        # Convert all lists stored under the 'template' key into tuples
        for key, list_ in v1['template'].items():
            v1['template'][key] = tuple(list_)
    return _recursive_mapping_proxy(templates)


#: An immutable mapping with additional values for ``CAT/data/workflow_dicts/workflow_yaml.yaml``.
#: Contains values which are generally not as easily represented in the .yaml format.
WORKFLOW_TEMPLATE: Mapping[str, _TemplateMapping] = finalize_templates()
