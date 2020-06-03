"""Tests for :mod:`CAT.data_handling.validation_schemas`."""

import os
from os.path import join

from schema import SchemaError
from unittest import mock

from scm.plams import AMSJob, ADFJob, Settings, CRSJob
from assertionlib import assertion

from CAT.utils import get_template
from CAT.data_handling.str_to_func import str_to_func
from CAT.data_handling.validation_schemas import (
    mol_schema, core_schema, ligand_schema, qd_schema, database_schema,
    mongodb_schema, bde_schema, qd_opt_schema, crs_schema, subset_schema
)

PATH = join('tests', 'test_files')


def test_mol_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.mol_schema`."""
    mol_dict = {}

    assertion.eq(mol_schema.validate(mol_dict), {'guess_bonds': False})

    mol_dict['guess_bonds'] = 1  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['guess_bonds'] = True  # Correct
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['is_core'] = 1  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['is_core'] = True
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['column'] = -1  # Exception: value < 0
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['column'] = 1.0  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['column'] = 1
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['row'] = -1  # Exception: value < 0
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['row'] = 1.0  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['row'] = 1
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['indices'] = 1.0  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['indices'] = [1, 5, 6, 7.0]  # Exception: an element has an incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['indices'] = (i for i in range(10))  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['indices'] = -1  # Exception: value < 0
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['indices'] = [-1, -2, -3, -4, -5]  # Exception: an element is < 0
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['indices'] = [1, 1, 2]  # Exception: duplicate elements
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)

    mol_dict['indices'] = 1
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1,))
    mol_dict['indices'] = [1, 2, 3, 4, 5]
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = {1, 2, 3, 4, 5}
    assertion.eq(mol_schema.validate(mol_dict)['indices'], (1, 2, 3, 4, 5))
    mol_dict['indices'] = (1, 2, 3, 4, 5)

    mol_dict['type'] = 1  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['type'] = 'bob'
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)

    mol_dict['name'] = 1  # Exception: incorrect type
    assertion.assert_(mol_schema.validate, mol_dict, exception=SchemaError)
    mol_dict['name'] = 'bob'
    assertion.eq(mol_schema.validate(mol_dict), mol_dict)


def test_database_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.database_schema`."""
    db_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'read': ('core', 'ligand', 'qd'),
        'write': ('core', 'ligand', 'qd'),
        'overwrite': (),
        'mongodb': {},
        'mol_format': ('pdb', 'xyz')

    }

    assertion.eq(database_schema.validate(db_dict), ref)

    for key in ('read', 'write', 'overwrite', 'mol_format', 'mongodb'):
        _db_dict = db_dict.copy()

        _db_dict[key] = 1  # Exception: incorrect type
        assertion.assert_(database_schema.validate, _db_dict, exception=SchemaError)
        _db_dict[key] = 'bob'  # Exception: incorrect value
        assertion.assert_(database_schema.validate, _db_dict, exception=SchemaError)
        _db_dict[key] = [1]  # Exception: element has incorrect type
        assertion.assert_(database_schema.validate, _db_dict, exception=SchemaError)
        _db_dict[key] = ['bob']  # Exception: element has incorrect value
        assertion.assert_(database_schema.validate, _db_dict, exception=SchemaError)

    db_dict['mongodb'] = True  # Exception: incorrect value
    assertion.assert_(database_schema.validate, db_dict, exception=SchemaError)
    db_dict['mongodb'] = False
    assertion.eq(database_schema.validate(db_dict), ref)
    db_dict['mongodb'] = {}
    assertion.eq(database_schema.validate(db_dict), ref)


def test_ligand_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.ligand_schema`."""
    lig_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'functional_groups': None,
        'optimize': {'job1': None},
        'split': True,
        'cosmo-rs': False,
        'cdft': False
    }

    assertion.eq(ligand_schema.validate(lig_dict), ref)

    lig_dict['optimize'] = 1  # Exception: incorrect type
    assertion.assert_(ligand_schema.validate, lig_dict, exception=SchemaError)
    lig_dict['optimize'] = True

    lig_dict['split'] = 1  # Exception: incorrect type
    assertion.assert_(ligand_schema.validate, lig_dict, exception=SchemaError)
    lig_dict['split'] = True

    lig_dict['cosmo-rs'] = 1  # Exception: incorrect type
    assertion.assert_(ligand_schema.validate, lig_dict, exception=SchemaError)
    lig_dict['cosmo-rs'] = {}
    assertion.eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {})
    lig_dict['cosmo-rs'] = False
    assertion.is_(ligand_schema.validate(lig_dict)['cosmo-rs'], False)
    lig_dict['cosmo-rs'] = True
    assertion.eq(ligand_schema.validate(lig_dict)['cosmo-rs'], {'job1': 'AMSJob'})

    lig_dict['functional_groups'] = 1  # Exception: incorrect type
    assertion.assert_(ligand_schema.validate, lig_dict, exception=SchemaError)
    lig_dict['functional_groups'] = 'CO'
    assertion.eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO']
    assertion.eq(ligand_schema.validate(lig_dict)['functional_groups'], ('CO',))
    lig_dict['functional_groups'] = ['CO', 'CO']  # Exception: duplicate elements
    assertion.assert_(ligand_schema.validate, lig_dict, exception=SchemaError)


def test_core_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.core_schema`."""
    core_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'dummy': 17,
        'allignment': 'sphere',
        'subset': None
    }

    assertion.eq(core_schema.validate(core_dict), ref)

    core_dict['dummy'] = 1.1  # Exception: incorrect value
    assertion.assert_(core_schema.validate, core_dict, exception=SchemaError)
    core_dict['dummy'] = 'H'
    assertion.eq(core_schema.validate(core_dict)['dummy'], 1)
    core_dict['dummy'] = 1
    assertion.eq(core_schema.validate(core_dict)['dummy'], 1)
    core_dict['dummy'] = 1.0
    assertion.eq(core_schema.validate(core_dict)['dummy'], 1)

    core_dict['allignment'] = 1.1  # Exception: incorrect type
    assertion.assert_(core_schema.validate, core_dict, exception=SchemaError)
    core_dict['allignment'] = 'bob'  # Exception: incorrect value
    assertion.assert_(core_schema.validate, core_dict, exception=SchemaError)
    core_dict['allignment'] = 'SPHERE'
    assertion.eq(core_schema.validate(core_dict)['allignment'], 'sphere')
    core_dict['allignment'] = 'surface'
    assertion.eq(core_schema.validate(core_dict)['allignment'], 'surface')


def test_qd_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_schema`."""
    qd_dict = {'dirname': '.'}
    ref = {
        'dirname': '.',
        'activation_strain': False,
        'optimize': False,
        'dissociate': False,
        'bulkiness': False,
        'construct_qd': True,
        'multi_ligand': None
    }

    assertion.eq(qd_schema.validate(qd_dict), ref)

    qd_dict['activation_strain'] = 1  # Exception: incorrect type
    assertion.assert_(qd_schema.validate, qd_dict, exception=SchemaError)
    qd_dict['activation_strain'] = True

    qd_dict['bulkiness'] = 1  # Exception: incorrect type
    assertion.assert_(qd_schema.validate, qd_dict, exception=SchemaError)
    qd_dict['bulkiness'] = False

    qd_dict['optimize'] = 1  # Exception: incorrect type
    assertion.assert_(qd_schema.validate, qd_dict, exception=SchemaError)
    qd_dict['optimize'] = True
    assertion.eq(qd_schema.validate(qd_dict)['optimize'], {'job1': 'AMSJob'})
    qd_dict['optimize'] = False
    assertion.is_(qd_schema.validate(qd_dict)['optimize'], False)

    qd_dict['dissociate'] = 1  # Exception: incorrect type
    assertion.assert_(qd_schema.validate, qd_dict, exception=SchemaError)
    qd_dict['dissociate'] = True  # Exception: incorrect value
    assertion.assert_(qd_schema.validate, qd_dict, exception=SchemaError)
    qd_dict['dissociate'] = False
    assertion.is_(qd_schema.validate(qd_dict)['dissociate'], False)


def test_mongodb_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.mongodb_schema`."""
    mongodb_dict = {}
    ref = {
        'host': 'localhost',
        'port': 27017
    }

    assertion.eq(mongodb_schema.validate(mongodb_dict), ref)

    mongodb_dict['port'] = 5.1  # Exception: incorrect value
    assertion.assert_(mongodb_schema.validate, mongodb_dict, exception=SchemaError)
    mongodb_dict['port'] = 5.0
    assertion.eq(mongodb_schema.validate(mongodb_dict)['port'], 5)
    mongodb_dict['port'] = 5
    assertion.eq(mongodb_schema.validate(mongodb_dict)['port'], 5)
    mongodb_dict['port'] = 27017

    mongodb_dict['host'] = 5.2  # Exception: incorrect type
    assertion.assert_(mongodb_schema.validate, mongodb_dict, exception=SchemaError)
    mongodb_dict['host'] = 'localhost'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['host'], 'localhost')
    mongodb_dict['host'] = 51
    assertion.eq(mongodb_schema.validate(mongodb_dict)['host'], 51)
    mongodb_dict['host'] = 51.0
    assertion.eq(mongodb_schema.validate(mongodb_dict)['host'], 51)

    mongodb_dict['username'] = 5.2  # Exception: incorrect type
    assertion.assert_(mongodb_schema.validate, mongodb_dict, exception=SchemaError)
    mongodb_dict['username'] = 'bob'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['username'], 'bob')
    mongodb_dict['username'] = 52
    assertion.eq(mongodb_schema.validate(mongodb_dict)['username'], 52)
    mongodb_dict['username'] = 52.0
    assertion.eq(mongodb_schema.validate(mongodb_dict)['username'], 52)

    mongodb_dict['password'] = 5.2  # Exception: incorrect type
    assertion.assert_(mongodb_schema.validate, mongodb_dict, exception=SchemaError)
    mongodb_dict['password'] = 'secret'
    assertion.eq(mongodb_schema.validate(mongodb_dict)['password'], 'secret')
    mongodb_dict['password'] = 53
    assertion.eq(mongodb_schema.validate(mongodb_dict)['password'], 53)
    mongodb_dict['password'] = 53.0
    assertion.eq(mongodb_schema.validate(mongodb_dict)['password'], 53)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_qd_opt_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.qd_opt_schema`."""
    _qd_opt_s1_default = get_template('qd.yaml')['UFF']
    _qd_opt_s2_default = _qd_opt_s1_default

    qd_opt_dict = Settings()
    ref = Settings({
        'job1': AMSJob,
        's1': _qd_opt_s1_default,
        'job2': AMSJob,
        's2': _qd_opt_s2_default,
        'keep_files': True,
        'use_ff': False
    })

    assertion.eq(qd_opt_schema.validate(qd_opt_dict), ref)

    for job in ('job1', 'job2'):
        qd_opt_dict[job] = 1  # Exception: incorrect type
        assertion.assert_(qd_opt_schema.validate, qd_opt_dict, exception=SchemaError)
        qd_opt_dict[job] = int  # Exception: incorrect value
        assertion.assert_(qd_opt_schema.validate, qd_opt_dict, exception=SchemaError)
        qd_opt_dict[job] = 'bob'  # Exception: incorrect value
        assertion.assert_(qd_opt_schema.validate, qd_opt_dict, exception=SchemaError)
        qd_opt_dict[job] = 'ADFJob'
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = 'ADFJOB'
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)
        qd_opt_dict[job] = ADFJob
        assertion.is_(qd_opt_schema.validate(qd_opt_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        qd_opt_dict[s] = 1  # Exception: incorrect type
        assertion.assert_(qd_opt_schema.validate, qd_opt_dict, exception=SchemaError)
        qd_opt_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)
        qd_opt_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(qd_opt_schema.validate(qd_opt_dict)[s], ref)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_crs_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.crs_schema`."""
    _crs_s1_default = get_template('qd.yaml')['COSMO-MOPAC']
    _crs_s2_default = get_template('qd.yaml')['COSMO-RS activity coefficient']
    _crs_s2_default.update(get_template('crs.yaml')['MOPAC PM6'])

    crs_dict = Settings()
    ref = Settings({
        'keep_files': True,
        'job1': AMSJob,
        's1': _crs_s1_default,
        'job2': CRSJob,
        's2': _crs_s2_default
    })

    assertion.eq(crs_schema.validate(crs_dict), ref)

    crs_dict['keep_files'] = 1  # Exception: incorrect type
    assertion.assert_(crs_schema.validate, crs_dict, exception=SchemaError)
    crs_dict['keep_files'] = False
    assertion.is_(crs_schema.validate(crs_dict)['keep_files'], False)
    crs_dict['keep_files'] = True

    for job in ('job1', 'job2'):
        crs_dict[job] = 1  # Exception: incorrect type
        assertion.assert_(crs_schema.validate, crs_dict, exception=SchemaError)
        crs_dict[job] = int  # Exception: incorrect value
        assertion.assert_(crs_schema.validate, crs_dict, exception=SchemaError)
        crs_dict[job] = 'bob'  # Exception: incorrect value
        assertion.assert_(crs_schema.validate, crs_dict, exception=SchemaError)
        crs_dict[job] = 'ADFJob'
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = 'ADFJOB'
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)
        crs_dict[job] = ADFJob
        assertion.is_(crs_schema.validate(crs_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        crs_dict[s] = 1  # Exception: incorrect type
        assertion.assert_(crs_schema.validate, crs_dict, exception=SchemaError)
        crs_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(crs_schema.validate(crs_dict)[s], ref)
        crs_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(crs_schema.validate(crs_dict)[s], ref)


@mock.patch.dict(os.environ,
                 {'ADFBIN': 'a', 'ADFHOME': '2019', 'ADFRESOURCES': 'b', 'SCMLICENSE': 'c'})
def test_bde_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.bde_schema`."""
    _bde_s1_default = get_template('qd.yaml')['MOPAC']
    # _bde_s2_default = get_template('qd.yaml')['UFF']

    bde_dict = Settings({'core_atom': 'Cd', 'lig_count': 2})
    ref = Settings({
        'use_ff': False,
        'core_atom': 48,
        'core_index': None,
        'core_core_dist': None,
        'lig_count': 2,
        'lig_core_pairs': 1,
        'lig_core_dist': None,
        'topology': None,

        'keep_files': True,
        'job1': AMSJob,
        's1': _bde_s1_default,
        'job2': None,
        's2': None
    })

    assertion.eq(bde_schema.validate(bde_dict), ref)

    bde_dict['keep_files'] = 1  # Exception: incorrect type
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['keep_files'] = False
    assertion.is_(bde_schema.validate(bde_dict)['keep_files'], False)
    bde_dict['keep_files'] = True

    bde_dict['core_atom'] = 5.1  # Exception: incorrect value
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['core_atom'] = 'H'
    assertion.eq(bde_schema.validate(bde_dict)['core_atom'], 1)
    bde_dict['core_atom'] = 1
    assertion.eq(bde_schema.validate(bde_dict)['core_atom'], 1)
    bde_dict['core_atom'] = 1.0
    assertion.eq(bde_schema.validate(bde_dict)['core_atom'], 1)

    bde_dict['lig_count'] = 5.2  # Exception: incorrect value
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['lig_count'] = -1  # Exception: incorrect value
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['lig_count'] = 3
    assertion.eq(bde_schema.validate(bde_dict)['lig_count'], 3)
    bde_dict['lig_count'] = 3.0
    assertion.eq(bde_schema.validate(bde_dict)['lig_count'], 3)

    bde_dict['core_index'] = 5.0  # Exception: incorrect type
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['core_index'] = [1, 2, 3, 4, 5.0]  # Exception: incorrect element type
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['core_index'] = [1, 2, 3, 4, 4]  # Exception: duplicate elements
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['core_index'] = 1
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1,))
    bde_dict['core_index'] = [1, 2, 3]
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))
    bde_dict['core_index'] = {1, 2, 3}
    assertion.eq(bde_schema.validate(bde_dict)['core_index'], (1, 2, 3))

    bde_dict['topology'] = 5.0  # Exception: incorrect type
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['topology'] = {'key': 'value'}  # Exception: incorrect value
    assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
    bde_dict['topology'] = {1: 'value'}
    assertion.eq(bde_schema.validate(bde_dict)['topology'], {1: 'value'})

    for dist in ('core_core_dist', 'lig_core_dist'):
        bde_dict[dist] = 'bob'  # Exception: incorrect type
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[dist] = -1  # Exception: incorrect value
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[dist] = 4
        assertion.eq(bde_schema.validate(bde_dict)[dist], 4.0)
        bde_dict[dist] = 4.0
        assertion.eq(bde_schema.validate(bde_dict)[dist], 4.0)

    for job in ('job1', 'job2'):
        bde_dict[job] = 1  # Exception: incorrect type
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[job] = int  # Exception: incorrect value
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[job] = 'bob'  # Exception: incorrect value
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[job] = 'ADFJob'
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = 'ADFJOB'
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)
        bde_dict[job] = ADFJob
        assertion.is_(bde_schema.validate(bde_dict)[job], ADFJob)

    ref = {'key1': {'key2': {'key3': True}}}

    for s in ('s1', 's2'):
        bde_dict[s] = 1  # Exception: incorrect type
        assertion.assert_(bde_schema.validate, bde_dict, exception=SchemaError)
        bde_dict[s] = {'key1': {'key2': {'key3': True}}}
        assertion.eq(bde_schema.validate(bde_dict)[s], ref)
        bde_dict[s] = join(PATH, 'settings.yaml')
        assertion.eq(bde_schema.validate(bde_dict)[s], ref)


def test_subset_schema() -> None:
    """Test :data:`CAT.data_handling.validation_schemas.subset_schema`."""
    subset_dict = {'f': 0.5}
    ref = {
        'f': 0.5,
        'mode': 'uniform',
        'start': None,
        'follow_edge': False,
        'cluster_size': 1,
        'randomness': None
    }

    subset_dict = subset_schema.validate(subset_dict)
    weight = subset_dict.pop('weight')

    assertion.eq(subset_dict, ref)
    assertion.function_eq(weight, str_to_func('np.exp(-x)'))

    subset_dict['p'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['p'] = 0  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['p'] = 0.33
    assertion.eq(subset_schema.validate(subset_dict)['p'], 0.33)
    subset_dict['p'] = 1
    assertion.eq(subset_schema.validate(subset_dict)['p'], 1)
    subset_dict['p'] = 999
    assertion.eq(subset_schema.validate(subset_dict)['p'], 999)
    subset_dict['p'] = -42.5
    assertion.eq(subset_schema.validate(subset_dict)['p'], -42.5)
    del subset_dict['p']

    subset_dict['weight'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['weight'] = 'np.exp(-y)'  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['weight'] = 'bob.exp(-x)'  # Exception: ImportError
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['weight'] = 'x.astype(int)'  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['weight'] = 'np.linalg.norm(x)'  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['weight'] = 'np.exp(-x)'
    subset_schema.validate(subset_dict)
    subset_dict['weight'] = 'x**2'
    subset_schema.validate(subset_dict)
    subset_dict['weight'] = 'scipy.special.expit(x)'
    subset_schema.validate(subset_dict)
    subset_dict['weight'] = '1 / (1 + np.exp(-x))'
    subset_schema.validate(subset_dict)
    subset_dict['weight'] = 'a = x**2; b = 5 * a; numpy.exp(b)'
    subset_schema.validate(subset_dict)
    del subset_dict['weight']

    subset_dict['f'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['f'] = -1  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['f'] = 0  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['f'] = 1.5  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['f'] = 0.33
    assertion.eq(subset_schema.validate(subset_dict)['f'], 0.33)
    subset_dict['f'] = 1
    assertion.eq(subset_schema.validate(subset_dict)['f'], 1.0)

    subset_dict['mode'] = 1  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['mode'] = 'bob'  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['mode'] = 'random'
    assertion.eq(subset_schema.validate(subset_dict)['mode'], 'random')
    subset_dict['mode'] = 'uniform'
    assertion.eq(subset_schema.validate(subset_dict)['mode'], 'uniform')
    subset_dict['mode'] = 'cluster'
    assertion.eq(subset_schema.validate(subset_dict)['mode'], 'cluster')

    subset_dict['start'] = 1.0  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['start'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['start'] = None
    assertion.is_(subset_schema.validate(subset_dict)['start'], None)
    subset_dict['start'] = 42
    assertion.eq(subset_schema.validate(subset_dict)['start'], 42)
    subset_dict['start'] = -5
    assertion.eq(subset_schema.validate(subset_dict)['start'], -5)
    subset_dict['start'] = 0
    assertion.eq(subset_schema.validate(subset_dict)['start'], 0)

    subset_dict['follow_edge'] = 1.0  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['follow_edge'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['follow_edge'] = True
    assertion.is_(subset_schema.validate(subset_dict)['follow_edge'], True)
    subset_dict['follow_edge'] = False
    assertion.is_(subset_schema.validate(subset_dict)['follow_edge'], False)

    subset_dict['cluster_size'] = 1.2  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = (i for i in range(10))  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = 0  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = (0,)  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = (None,)  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['cluster_size'] = 10
    assertion.eq(subset_schema.validate(subset_dict)['cluster_size'], 10)
    subset_dict['cluster_size'] = 10.0
    assertion.eq(subset_schema.validate(subset_dict)['cluster_size'], 10)
    subset_dict['cluster_size'] = [1, 5, 10]
    assertion.eq(subset_schema.validate(subset_dict)['cluster_size'], (1, 5, 10))
    subset_dict['cluster_size'] = [1.0, 5.0, 10]
    assertion.eq(subset_schema.validate(subset_dict)['cluster_size'], (1, 5, 10))

    subset_dict['randomness'] = 'bob'  # Exception: incorrect type
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['randomness'] = -1  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['randomness'] = 1.5  # Exception: incorrect value
    assertion.assert_(subset_schema.validate, subset_dict, exception=SchemaError)
    subset_dict['randomness'] = 0.33
    assertion.eq(subset_schema.validate(subset_dict)['randomness'], 0.33)
    subset_dict['randomness'] = 1
    assertion.eq(subset_schema.validate(subset_dict)['randomness'], 1.0)
    subset_dict['randomness'] = 0
    assertion.eq(subset_schema.validate(subset_dict)['randomness'], 0)
