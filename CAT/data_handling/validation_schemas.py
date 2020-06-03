"""
CAT.data_handling.validation_schemas
====================================

A module designed for sanitizing and interpreting the input file.

Index
-----
.. currentmodule:: CAT.data_handling.validation_schemas
.. autosummary::
    mol_schema
    core_schema
    ligand_schema
    qd_schema
    database_schema
    mongodb_schema
    ligand_opt_schema
    bde_schema
    qd_opt_schema
    crs_schema
    asa_schema
    subset_schema
    multi_ligand_schema
    JOB_MAP
    JOB_MAP_SCM

API
---
.. autodata:: mol_schema
    :annotation: : schema.Schema
.. autodata:: core_schema
    :annotation: : schema.Schema
.. autodata:: ligand_schema
    :annotation: : schema.Schema
.. autodata:: qd_schema
    :annotation: : schema.Schema
.. autodata:: database_schema
    :annotation: : schema.Schema
.. autodata:: mongodb_schema
    :annotation: : schema.Schema
.. autodata:: ligand_opt_schema
    :annotation: : schema.Schema
.. autodata:: bde_schema
    :annotation: : schema.Schema
.. autodata:: qd_opt_schema
    :annotation: : schema.Schema
.. autodata:: asa_schema
    :annotation: : schema.Schema
.. autodata:: crs_schema
    :annotation: : schema.Schema
.. autodata:: subset_schema
    :annotation: : schema.Schema
.. autodata:: multi_ligand_schema
    :annotation: : schema.Schema
.. autodata:: JOB_MAP
    :annotation: : Mapping[str, Type[Job]]
.. autodata:: JOB_MAP_SCM
    :annotation: : Mapping[str, Type[Job]]

"""

from types import MappingProxyType
from typing import Collection, Callable, Any, Optional, TypeVar, Mapping, Type, Tuple
from operator import __index__
from collections import abc

import numpy as np
from schema import Or, And, Use, Schema
from schema import Optional as Optional_

from scm.plams import CRSJob, Settings
from scm.plams.core.basejob import Job

from scm.plams.interfaces.adfsuite.adf import ADFJob
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.interfaces.adfsuite.uff import UFFJob
from scm.plams.interfaces.adfsuite.band import BANDJob
from scm.plams.interfaces.adfsuite.dftb import DFTBJob
from scm.plams.interfaces.adfsuite.mopac import MOPACJob
from scm.plams.interfaces.adfsuite.reaxff import ReaxFFJob

from scm.plams.interfaces.thirdparty.cp2k import Cp2kJob
from scm.plams.interfaces.thirdparty.orca import ORCAJob
from scm.plams.interfaces.thirdparty.dirac import DiracJob
from scm.plams.interfaces.thirdparty.gamess import GamessJob
from scm.plams.interfaces.thirdparty.dftbplus import DFTBPlusJob

from .str_to_func import str_to_func
from ..utils import get_template, validate_path, validate_core_atom, check_sys_var
from ..mol_utils import to_atnum

__all__ = ['mol_schema', 'core_schema', 'ligand_schema', 'qd_schema', 'database_schema',
           'mongodb_schema', 'bde_schema', 'ligand_opt_schema', 'qd_opt_schema', 'crs_schema',
           'asa_schema', 'subset_schema']


def val_float(value: Any) -> bool:
    """Check if a float-like object has been passed (:data:`typing.SupportsFloat`)."""
    try:
        value.__float__()
        return True
    except Exception:
        return False


def val_int(value: Any) -> bool:
    """Check if a int-like object has been passed (:data:`typing.SupportsInt`)."""
    try:
        value.__int__()
        return float(value).is_integer()
    except Exception:
        return False


def val_index(value: Any) -> bool:
    """Check if a index-like object has been passed (:data:`typing.SupportsIndex`)."""
    try:
        value.__index__()
        return True
    except Exception:
        return False


#: A dictionary for translating strings into :class:`plams.Job` types for third party software
JOB_MAP: Mapping[str, Type[Job]] = MappingProxyType({
    'cp2k': Cp2kJob, 'cp2kjob': Cp2kJob,
    'orca': ORCAJob, 'orcajob': ORCAJob,
    'dirac': DiracJob, 'diracjob': DiracJob,
    'gamess': GamessJob, 'gamessjob': GamessJob,
    'dftbplus': DFTBPlusJob, 'dftbplusjob': DFTBPlusJob,
})


#: A dictionary for translating strings into :class:`plams.Job` types for ADF
JOB_MAP_SCM: Mapping[str, Type[Job]] = MappingProxyType({
    'adf': ADFJob, 'adfjob': ADFJob,
    'ams': AMSJob, 'amsjob': AMSJob,
    'uff': UFFJob, 'uffjob': UFFJob,
    'band': BANDJob, 'bandjob': BANDJob,
    'dftb': DFTBJob, 'dftbjob': DFTBJob,
    'mopac': MOPACJob, 'mopacjob': MOPACJob,
    'reaxff': ReaxFFJob, 'reaxffjob': ReaxFFJob,
    'crs': CRSJob, 'cosmo-rs': CRSJob, 'crsjob': CRSJob
})


DB_NAMES: Tuple[str, ...] = ('core', 'ligand', 'qd')
FORMAT_NAMES: Tuple[str, ...] = ('pdb', 'xyz')
FORMAT_NAMES2: Tuple[str, ...] = ('pdb', 'xyz', 'mol', 'mol2')


def val_job_type(value: type) -> type:
    """Call :func:`.check_sys_var` if value is in :data:`.JOB_MAP_SCM`."""
    if value in JOB_MAP_SCM.values():
        check_sys_var()
    return value


def str_to_job_type(key: str) -> type:
    """Convert a string into a type object.

    Parameters
    ----------
    key : str
        An alias for a :class:`type` object of a PLAMS :class:`Job`.
        **key** is converted by passing it into :data:`.JOB_MAP` or :data:`.JOB_MAP_SCM`.
        The :func:`check_sys_var` function is called if **key** is present
        in :data:`.JOB_MAP_SCM`.

    Returns
    -------
    |type|_
        A type object constructed from **key**.

    Raises
    ------
    KeyError
        Raised if the decapitalized **key** is available in neither :data:`.JOB_MAP` nor
        :data:`.JOB_MAP_SCM`.

    """
    _key = key.lower()
    if _key in JOB_MAP:
        return JOB_MAP[_key]
    elif _key in JOB_MAP_SCM:
        check_sys_var()
        return JOB_MAP_SCM[_key]
    raise KeyError(f'No Job type alias available for {repr(_key)}')


T = TypeVar('T')


def to_tuple(collection: Collection[T], func: Optional[Callable[[T], Any]] = None) -> tuple:
    """Convert a collection to a sorted tuple."""
    try:
        ret = sorted(collection)
    except TypeError:  # The collection contains a mix of sorting-incompatible objects
        ret = sorted(collection, key=str)
    finally:
        if func is None:
            return tuple(ret)
        else:
            return tuple(func(i) for i in ret)


# The **default** parameter of schema.Optional() will automatically call any callable
# Solution: provide a callable that returns another callable
def _get_amsjob() -> type:
    """Return a type object of :class:`.AMSJob`."""
    return AMSJob


def _get_crsjob() -> type:
    """Return a type object of :class:`.CRSJob`."""
    return CRSJob


# Default settings templates
_bde_s1_default = get_template('qd.yaml')['MOPAC']
_bde_s2_default = get_template('qd.yaml')['UFF']

_qd_opt_s1_default = get_template('qd.yaml')['UFF']
_qd_opt_s2_default = _qd_opt_s1_default

_crs_s1_default = get_template('qd.yaml')['COSMO-MOPAC']
_crs_s2_default = get_template('qd.yaml')['COSMO-RS activity coefficient']
_crs_s2_default.update(get_template('crs.yaml')['MOPAC PM6'])


#: Schema for validating the ``['input_ligands']`` and ``['input_cores']`` blocks.
mol_schema: Schema = Schema({
    Optional_('guess_bonds', default=False):
        And(bool, error=".guess_bonds expects a boolean"),

    Optional_('is_core'):
        And(bool, error=".is_core expects a boolean"),

    Optional_('is_qd'):
        And(bool, error=".is_qd expects a boolean"),

    Optional_('ligand_smiles'):
        And(str, error=".ligand_smiles expects a string"),

    Optional_('ligand_anchor'):
        And(str, error=".ligand_smiles expects a string"),

    Optional_('column'):
        And(val_index, lambda n: n.__index__() >= 0, Use(int),
            error=".column expects an integer larger than or equal to 0"),

    Optional_('row'):
        And(val_index, lambda n: n.__index__() >= 0, Use(int),
            error=".row expects an integer larger than or equal to 0"),

    Optional_('indices'):
        Or(
            And(
                val_index, lambda n: n.__index__() >= 0,
                Use(lambda n: to_tuple([n], func=__index__)),
                error=".indices expects an integer larger than or equal to 0"
            ),
            And(
                abc.Collection,
                lambda n: all(val_index(i) for i in n),
                lambda n: all(i.__index__() >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(lambda n: to_tuple(n, func=__index__)),
                error=".indices expects one or more unique integers larger than or equal to 0"
            ),
            error=".indices expects an atomic index (int) or a list unique atomic indices"
        ),

    Optional_('type'):
        And(str, error='.type expects a string'),

    Optional_('name'):
        And(str, error='.name expects a string'),

    Optional_('path'):
        And(
            Use(validate_path),
            error=".path expects a None or a string pointing to an existing directory"
        ),
})

#: Schema for validating the ``['optional']['core']`` block.
core_schema: Schema = Schema({
    'dirname':
        And(str, error='optional.core.dirname expects a string'),

    Optional_('dummy', default=17):  # Return a tuple of atomic numbers
        Or(
            And(val_int, Use(lambda n: to_atnum(int(n)))),
            And(str, Use(to_atnum)),
            error='optional.core.dummy expects a valid atomic number (int) or symbol (string)'
        ),

    Optional_('subset', default=None):
        Or(None, dict, error="optional.core.subset epected 'None' or a dictionary"),

    Optional_('allignment', default='sphere'):
        And(str, lambda n: n.lower() in {'sphere', 'surface'},
            Use(str.lower), error="optional.core.allignment expected 'sphere' or 'surface'")
})

#: Schema for validating the ``['optional']['core']['subset']`` block.
subset_schema: Schema = Schema({
    'f':
        And(val_float, lambda n: 0 < float(n) <= 1, Use(float)),

    Optional_('p'):
        And(val_float, lambda n: abs(float(n)) > 0, Use(float)),

    Optional_('mode', default='uniform'):
        And(str, lambda n: n.lower() in {'uniform', 'random', 'cluster'}, Use(str.lower)),

    Optional_('start', default=None):
        Or(
            None,
            And(val_index, Use(__index__))
        ),

    Optional_('follow_edge', default=False):
        bool,

    Optional_('cluster_size', default=1):
        Or(
            And(val_int, lambda n: int(n) > 0, Use(int)),
            And(abc.Collection,
                lambda n: all(val_int(i) for i in n),
                lambda n: all(int(i) > 0 for i in n),
                Use(lambda n: to_tuple(n, func=int)))
        ),

    Optional_('weight', default=lambda: str_to_func('np.exp(-x)')):
        Or(
            And(Callable, Use(str_to_func)),
            And(str, lambda n: 'x' in n, Use(str_to_func))
        ),

    Optional_('randomness', default=None):
        Or(
            None,
            And(val_float, lambda n: 0 <= float(n) <= 1, Use(float))
        )

})

#: Schema for validating the ``['optional']['database']`` block.
database_schema: Schema = Schema({
    # path+directory name of the database
    'dirname':
        And(str, error='optional.database.dirname expects a string'),

    Optional_('read', default=DB_NAMES):  # Attempt to pull structures from the database
        Or(
            And(None, Use(lambda n: ())),
            And(bool, Use(lambda n: DB_NAMES if n is True else ())),
            And(
                str,
                lambda n: n in DB_NAMES,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.read are: {repr(DB_NAMES)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in DB_NAMES for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.read are: {repr(DB_NAMES)}'
            ),
            error='optional.database.read expects a boolean, string or list of unique strings'
        ),

    Optional_('write', default=DB_NAMES):  # Attempt to write structures to the database
        Or(
            And(None, Use(lambda n: ())),
            And(bool, Use(lambda n: DB_NAMES if n is True else ())),
            And(
                str,
                lambda n: n in DB_NAMES,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.write are: {repr(DB_NAMES)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in DB_NAMES for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.write are: {repr(DB_NAMES)}'
            ),
            error='optional.database.write expects a boolean, string or list of unique strings'
        ),

    Optional_('overwrite', default=tuple):  # Allow previous entries to be overwritten
        Or(
            And(None, Use(lambda n: ())),
            And(bool, Use(lambda n: DB_NAMES if n is True else ())),
            And(
                str,
                lambda n: n in DB_NAMES,
                Use(lambda n: (n,)),
                error=f'allowed values for optional.database.overwrite are: {repr(DB_NAMES)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in DB_NAMES for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.overwrite are: {repr(DB_NAMES)}'
            ),
            error='optional.database.overwrite expects a boolean, string or list of unique strings'
        ),

    Optional_('mongodb', default=dict):  # Settings specific to MongoDB
        Or(
            dict,
            And(bool, lambda n: n is False, Use(lambda n: {})),
            error='optional.database.mongodb expects True (boolean) or a dictionary'
        ),

    Optional_('mol_format', default=FORMAT_NAMES):  # Return a tuple of file formats
        Or(
            And(None, Use(lambda n: ())),
            And(bool, Use(lambda n: FORMAT_NAMES if n is True else ())),
            And(
                str,
                lambda n: n in FORMAT_NAMES,
                error=f'allowed values for optional.database.mol_format are: {repr(FORMAT_NAMES2)}'
            ),
            And(
                abc.Collection,
                lambda n: all(i in FORMAT_NAMES for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error=f'allowed values for optional.database.mol_format are: {repr(FORMAT_NAMES2)}'
            ),
            error='optional.database.mol_format expects a boolean, string or list of unique strings'
        )
})


#: Schema for validating the ``['optional']['ligand']`` block.
ligand_schema: Schema = Schema({
    # path+directory name of the ligand directory
    'dirname':
        And(str, error='optional.ligand.dirname expects a string'),

    Optional_('functional_groups', default=None):
        Or(
            None,
            And(str, Use(lambda n: (n,))),
            And(
                abc.Collection,
                lambda n: all(isinstance(i, str) for i in n),
                lambda n: len(n) == len(set(n)),
                Use(to_tuple),
                error='optional.ligand.functional_groups expects a list of unique SMILES strings'
            ),
            error=('optional.ligand.functional_groups expects None (NoneType), a SMILES string, '
                   'or a list of unique SMILES string')
        ),

    Optional_('optimize', default={'job1': None}):  # Optimize the ligands
        Or(
            And(
                bool, Use(lambda n: ({'job1': None} if n else False)),
                error='optional.ligand.cosmo-rs expects a boolean or dictionary'
            ),
            And(
                dict,
                error='optional.ligand.cosmo-rs expects a boolean or dictionary'
            )
        ),

    Optional_('split', default=True):  # Remove a counterion from the function group
        And(bool, error='optional.ligand.split expects a boolean'),

    Optional_('cosmo-rs', default=False):  # Settings specific to ligand COSMO-RS calculations
        Or(
            dict,
            And(bool, Use(lambda n: {'job1': 'AMSJob'} if n else False)),
            error='optional.ligand.cosmo-rs expects a boolean or dictionary'
        ),

    Optional_('cdft', default=False):  # Settings specific to ligand conceptual dft calculations
        Or(
            dict,
            And(bool, Use(lambda n: {'job1': 'ADFJob'} if n else False)),
            error='optional.ligand.cdft expects a boolean or dictionary'
        ),
})


#: Schema for validating the ``['optional']['qd']`` block.
qd_schema: Schema = Schema({
    # path+directory name of the quantum dot directory
    'dirname':
        And(str, error='optional.qd.dirname expects a string'),

    Optional_('construct_qd', default=True):  # Construct quantum dots
        Or(
            bool,
            error='optional.qd.construct_qd expects a boolean'
        ),

    # Settings specific to a quantum dot activation strain analyses
    Optional_('activation_strain', default=False):
        Or(
            And(bool, Use(lambda n: {'job1': None} if n else False)),
            dict, error='optional.qd.activation_strain expects a boolean or dictionary'
        ),

    Optional_('bulkiness', default=False):  # Ligand bulkiness workflow
        And(bool, error='optional.qd.bulkiness expects a boolean'),

    Optional_('optimize', default=False):  # Settings for quantum dot geometry optimizations
        Or(
            dict,
            And(bool, Use(lambda n: ({'job1': 'AMSJob'} if n else False))),
            error='optional.qd.optimize expects a boolean or dictionary'
        ),

    # Settings for quantum dot ligand dissociation calculations
    Optional_('dissociate', default=False):
        Or(
            dict,
            And(bool, lambda n: n is False),
            error='optional.qd.dissociate expects False (boolean) or a dictionary'
        ),

    Optional_('multi_ligand', default=None):  # Settings for multi-ligand attachment
        Or(
            None, dict,
            error='optional.qd.multi_ligand expects None or dictionary'
        ),
})


#: Schema for validating the ``['optional']['database']['mongodb']`` block.
mongodb_schema: Schema = Schema({
    # Optional username for the MongoDB host
    Optional_('username'):
        Or(str, And(val_int, Use(int)),
           error='optional.database.mongodb.username expects a string or integer'),

    Optional_('password'):  # Optional password for the MongoDB host
        Or(str, And(val_int, Use(int)),
           error='optional.database.mongodb.password expects a string or integer'),

    Optional_('host', default='localhost'):  # Name of the MongoDB host
        Or(str, And(val_int, Use(int)),
           error='optional.database.mongodb.host expects a string or integer'),

    Optional_('port', default=27017):  # Port of the MongoDB host
        And(val_int, Use(int), error='optional.database.mongodb.port expects an integer'),

    Optional_(str):  # Other keyword arguments for :class:`pymongo.MongoClient`
        object
})


#: Schema for validating the ``['optional']['ligand']['optimize']`` block.
ligand_opt_schema: Schema = Schema({
    Optional_('use_ff', default=False):
        bool,

    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        And(bool, error='optional.ligand.optimize.keep_files expects a boolean'),

    # The Job type and settings for the conformation search
    Optional_('job1', default=None): None,
    Optional_('s1', default=None): dict,

    # The Job type for the final geometry optimization
    Optional_('job2', default=None):
        Or(
            None,
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.optimize.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error=('optional.ligand.optimize.job1 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.ligand.optimize.job2 expects a string or a type object'
        ),

    # The Job Settings for the final geometry optimization
    Optional_('s2', default=Settings):
        Or(
            None,
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.optimize.s2 expects a string or a dictionary'
        ),
})


#: Schema for validating the ``['optional']['qd']['dissociate']`` block.
bde_schema: Schema = Schema({
    Optional_('use_ff', default=False):
        bool,

    # Atom type of the to-be dissociated core atom
    'core_atom':
        And(
            Or(And(val_int, Use(int)), str), Use(validate_core_atom),
            error=('optional.qd.dissociate.core_atom expects a SMILES string, '
                   'atomic number (int) or atomic symbol (str)')
        ),

    'lig_count':  # The number of ligands per core_atom
        And(
            val_int, lambda n: int(n) >= 0, Use(int),
            error='optional.qd.dissociate.lig_count expects an integer larger than or equal to 0'
        ),

    Optional_('keep_files', default=True):  # Delete files after the calculations are finished
        And(bool, error='optional.qd.dissociate.keep_files expects a boolean'),

    Optional_('core_core_dist', default=None):
        Or(
            None,
            And(
                val_float, lambda n: float(n) > 0.0, Use(float),
                error=('optional.qd.dissociate.core_core_dist expects an integer or float '
                       'larger than or equal to 0.0')
            )
        ),

    Optional_('lig_core_dist', default=None):
        Or(
            None,
            And(
                val_float, lambda n: float(n) > 0.0, Use(float),
                error=('optional.qd.dissociate.lig_core_dist expects an integer or float '
                       'larger than or equal to 0.0')
            )
        ),

    Optional_('lig_core_pairs', default=1):
        Or(
            None,
            And(
                val_int, lambda n: int(n) > 0, Use(int),
                error=('optional.qd.dissociate.lig_pairs expects an integer'
                       'larger than or equal to 0')
            )
        ),

    Optional_('core_index', default=None):
        Or(
            None,
            And(
                val_index, lambda n: n.__index__() >= 0,
                Use(lambda n: to_tuple([n], func=__index__)),
                error=('optional.qd.dissociate.core_index expects an integer '
                       'larger than or equal to 0')
            ),
            And(
                abc.Collection,
                lambda n: all(val_index(i) for i in n),
                lambda n: all(i.__index__() >= 0 for i in n),
                lambda n: len(n) == len(set(n)),
                Use(lambda n: to_tuple(n, func=__index__)),
                error=('optional.qd.dissociate.core_index expects a list of unique integers '
                       'larger than or equal to 0')
            ),
            error=('optional.qd.dissociate.core_index expects an integer or list of unique integers'
                   'larger than or equal to 0')
        ),

    Optional_('topology', default=None):
        Or(
            None,
            And(
                abc.Mapping, lambda n: all(val_int(k) for k in n.keys()),
                Use(lambda n: {int(k): v for k, v in n.items()}),
                error='optional.qd.dissociate.topology expects a dictionary with integers as keys'
                )
        ),


    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.dissociate.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error='optional.qd.dissociate.job1 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.dissociate.job1 expects a string or a type object'

        ),

    Optional_('s1', default=_bde_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.dissociate.s1 expects a string or a dictionary'
        ),

    Optional_('job2', default=None):
        Or(
            None,
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.dissociate.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error='optional.qd.dissociate.job2 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.dissociate.job2 expects a string or a type object'
        ),

    Optional_('s2', default=None):
        Or(
            None,
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.dissociate.s2 expects a string or a dictionary'
        )
})

#: Schema for validating the ``['optional']['qd']['optimize']`` block.
qd_opt_schema: Schema = Schema({
    Optional_('use_ff', default=False):
        bool,

    Optional_('keep_files', default=True):
        And(bool, error='optional.qd.opt.keep_files expects a boolean'),

    # The job type for the first half of the optimization
    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.opt.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error='optional.qd.opt.job1 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.opt.job1 expects a string or a type object'
        ),

    # The job settings for the first half of the optimization
    Optional_('s1', default=_qd_opt_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.opt.s1 expects a string or a dictionary'
        ),

    # The job type for the second half of the optimization
    Optional_('job2', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.opt.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error='optional.qd.opt.job2 expects a string that is a valid plams.Job alias'
            ),
            error='optional.qd.opt.job2 expects a string or a type object'
        ),

    # The job settings for the second half of the optimization
    Optional_('s2', default=_qd_opt_s2_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.opt.s2 expects a string or a dictionary'
        )
})

#: Schema for validating the ``['optional']['ligand']['cosmo-rs']`` block.
crs_schema: Schema = Schema({
    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        And(bool, error='optional.ligand.cosmo-rs.keep_files expects a boolean'),

    # The job type for constructing the COSMO surface
    Optional_('job1', default=_get_amsjob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.cosmo-rs.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error=('optional.ligand.cosmo-rs.job1 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.ligand.cosmo-rs.job1 expects a string or a type object'
        ),

    # The settings for constructing the COSMO surface
    Optional_('s1', default=_crs_s1_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.cosmo-rs.s1 expects a string or a dictionary'
        ),

    Optional_('job2', default=_get_crsjob):  # The job type for the actual COSMO-RS calculation
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.cosmo-rs.job2 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error=('optional.ligand.cosmo-rs.job2 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.ligand.cosmo-rs.job2 expects a string or a type object'
        ),

    # The settings for the actual COSMO-RS calculation
    Optional_('s2', default=_crs_s2_default.copy()):
        Or(
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.cosmo-rs.s2 expects a string or a dictionary'
        )
})

#: Schema for validating the ``['optional']['qd']['activation_strain']`` block.
asa_schema: Schema = Schema({
    Optional_('use_ff', default=False):
        bool,

    Optional_('md', default=False):
        bool,

    Optional_('dump_csv', default=False):
        bool,

    Optional_('iter_start', default=500):
        And(val_index, lambda n: n.__index__() >= 0, Use(__index__)),

    Optional_('el_scale14', default=1.0):
        And(val_float, Use(float)),

    Optional_('lj_scale14', default=1.0):
        And(val_float, Use(float)),

    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        And(bool, error='optional.qd.activation_strain.keep_files expects a boolean'),

    Optional_('distance_upper_bound', default=np.inf):
        Or(
            And(str, lambda n: 'inf' in n.lower(), Use(lambda n: np.inf)),
            And(val_float, lambda n: float(n) > 0, Use(float))
        ),

    Optional_('shift_cutoff', default=True):
        bool,

    Optional_('k', default=20):
        And(val_int, lambda n: int(n) > 0, Use(int)),

    # The job type for constructing the COSMO surface
    Optional_('job1', default=None):
        Or(
            None,
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.qd.activation_strain.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error=('optional.qd.activation_strain.job1 expects a string '
                       'that is a valid plams.Job alias')
            ),
            error='optional.qd.activation_strain.job1 expects a string or a type object'
        ),

    # The settings for constructing the COSMO surface
    Optional_('s1', default=None):
        Or(
            None,
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.qd.activation_strain.s1 expects a string or a dictionary'
        ),
})

#: Schema for validating the ``['optional']['qd']['multi_ligand']`` block.
multi_ligand_schema: Schema = Schema({
    'ligands':
        And(
            abc.Collection,
            lambda n: not isinstance(n, str),
            lambda n: all(isinstance(i, str) for i in n),
            Use(tuple)
        ),

    Optional_('dummy', default=None):
        Or(
            None,
            And(abc.Collection,
                lambda n: not isinstance(n, str),
                lambda n: len(set(n)) == len(n),
                Use(lambda n: to_tuple(n, func=to_atnum)))
        ),

    Optional_('f', default=None):
        Or(
            None,
            And(
                abc.Collection,
                lambda n: all(val_float(i) for i in n),
                lambda n: all(float(i) > 0 for i in n),
                lambda n: 0 < sum(float(i) for i in n) <= 1,
                Use(lambda n: to_tuple(n, func=float))
            )
        ),

    Optional_('mode', default='uniform'):
        And(str, lambda n: n.lower() in {'uniform', 'random', 'cluster'}, Use(str.lower)),

    Optional_('start', default=None):
        Or(
            None,
            And(val_index, Use(__index__))
        ),

    Optional_('follow_edge', default=False):
        bool,

    Optional_('cluster_size', default=1):
        Or(
            And(val_int, lambda n: int(n) > 0, Use(int)),
            And(abc.Collection,
                lambda n: all(val_int(i) for i in n),
                lambda n: all(int(i) > 0 for i in n),
                Use(lambda n: to_tuple(n, func=int)))
        ),

    Optional_('weight', default=lambda: str_to_func('np.exp(-x)')):
        Or(
            And(Callable, Use(str_to_func)),
            And(str, lambda n: 'x' in n, Use(str_to_func))
        ),

    Optional_('randomness', default=None):
        Or(
            None,
            And(val_float, lambda n: 0 <= float(n) <= 1, Use(float))
        )
})


#: Schema for validating the ``['optional']['ligand']['cdft']`` block.
cdft_schema: Schema = Schema({
    # Delete files after the calculations are finished
    Optional_('keep_files', default=True):
        And(bool, error='optional.ligand.cdft.keep_files expects a boolean'),

    # The Job type for the final geometry optimization
    Optional_('job1', default=lambda n: ADFJob):
        Or(
            And(
                And(type, lambda n: issubclass(n, Job), Use(val_job_type)),
                error=('optional.ligand.cdft.job1 expects a type object '
                       'that is a subclass of plams.Job')
            ),
            And(
                str, Use(str_to_job_type),
                error=('optional.ligand.cdft.job1 expects a string '
                       'that is a valid plams.Job alias')
            ),
        ),

    # The Job Settings for the final geometry optimization
    Optional_('s2', default=Settings):
        Or(
            None,
            dict,
            And(str, Use(lambda n: get_template(n, from_cat_data=False))),
            error='optional.ligand.cdft.s1 expects a string or a dictionary'
        ),
})
