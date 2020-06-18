"""Tests for :mod:`CAT.data_handling.entry_points`."""

import os
import warnings
from pathlib import Path

import numpy as np

from scm.plams import Molecule
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.data_handling.warn_map import MoleculeWarning
from CAT.data_handling.entry_points import main

PATH = Path('tests') / 'test_files'

LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'
PATH_REF = PATH / 'ligand_ref'


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_main() -> None:
    """Test :func:`CAT.data_handling.entry_points.main`."""
    filename = str(PATH / 'input2.yaml')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', MoleculeWarning)
        main([filename])

    assertion.assert_(main, [f'{filename}bob'], exception=FileNotFoundError)

    # Compare ligand geometries
    iterator = ((f, Molecule(LIG_PATH / f), Molecule(PATH_REF / f)) for
                f in os.listdir(PATH_REF) if f.endswith('xyz'))

    for f, mol, mol_ref in iterator:
        if not f.endswith('xyz'):
            continue
        np.testing.assert_allclose(mol, mol_ref, err_msg=f, rtol=0, atol=10**-2)
