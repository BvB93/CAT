"""Tests for :mod:`CAT.data_handling.mol_to_file`."""

from os import mkdir
from os.path import (join, isdir)
from shutil import rmtree

from scm.plams import Molecule
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.data_handling.mol_to_file import mol_to_file

_PATH = join('tests', 'test_files')
PATH = join(_PATH, 'mol_to_file')
MOL = Molecule(join(_PATH, 'Methanol.pdb'))


@delete_finally(PATH)
def test_mol_to_file() -> None:
    """Test :func:`CAT.data_handling.mol_to_file.mol_to_file`."""
    if not isdir(PATH):
        mkdir(PATH)

    _ref = ('mol.xyz', 'mol.pdb', 'mol.mol', 'mol.mol2')
    ref = [join(PATH, item) for item in _ref]

    mol = MOL.copy()
    mol.properties.name = 'mol'

    mol_list = [mol]
    mol_format = ('xyz', 'pdb', 'mol', 'mol2')
    kwargs = {'path': PATH, 'mol_format': mol_format}

    mol_to_file(mol_list, **kwargs)
    for file in ref:
        assertion.isfile(file)

    rmtree(PATH)
    mkdir(PATH)
    kwargs['mol_format'] = kwargs['mol_format'][0:2]
    mol_to_file(mol_list, **kwargs)
    for file in ref[:2]:
        assertion.isfile(file)
    for file in ref[2:]:
        assertion.isfile(file, invert=True)

    kwargs['path'] = join(PATH, 'bob')
    assertion.assert_(mol_to_file, mol_list, exception=FileNotFoundError, **kwargs)

    kwargs['path'] = join(PATH, 'mol.xyz')
    assertion.assert_(mol_to_file, mol_list, exception=NotADirectoryError, **kwargs)

    kwargs['path'] = PATH
    kwargs['mol_format'] = ('bob')
    assertion.assert_(mol_to_file, mol_list, exception=ValueError, **kwargs)
