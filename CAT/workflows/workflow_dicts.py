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
    V_BULK
)

if TYPE_CHECKING:
    if sys.version_info >= (3, 8):
        from typing import TypedDict
    else:
        from typing_extensions import TypedDict

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
