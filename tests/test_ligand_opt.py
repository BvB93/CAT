"""Tests for :mod:`CAT.attachment.ligand_opt`."""

from os.path import join

import numpy as np

from scm.plams import readpdb, to_rdmol
from CAT.attachment.ligand_opt import rdmol_as_array

PATH = join('tests', 'test_files')
MOL = readpdb(join(PATH, 'Methanol.pdb'))


def test_rdmol_as_array() -> None:
    """Test :func:`CAT.attachment.ligand_opt.rdmol_as_array`."""
    xyz_ref = np.array(MOL)

    rdmol = to_rdmol(MOL)
    xyz = rdmol_as_array(rdmol)
    np.testing.assert_allclose(xyz, xyz_ref)
