"""Tests for :mod:`CAT.mol_utils`."""

from os.path import join
from itertools import chain

from scm.plams import (Molecule, PeriodicTable, PTError)
import scm.plams.interfaces.molecule.rdkit as molkit

from CAT.assertion.assertion_manager import assertion
from CAT.mol_utils import (
    from_mol_other, from_rdmol, get_index, merge_mol, separate_mod,
    to_atnum, to_symbol, adf_connectivity, fix_carboxyl, fix_h
)

PATH = join('tests', 'test_files')
MOL = Molecule(join(PATH, 'Methanol.xyz'))  # Methanol; BP86/QZ4P
MOL.guess_bonds()


def test_from_mol_other() -> None:
    """Test :meth:`Molecule.from_mol_other`."""
    mol = MOL.copy()
    mol_rot = Molecule(join(PATH, 'Methanol_rotate.xyz'))
    mol.from_mol_other(mol_rot)

    for at, at_ref in zip(mol, mol_rot):
        assertion.eq(at.coords, at_ref.coords)
        assertion.eq(at.symbol, at_ref.symbol)


def test_from_rdmol() -> None:
    """Test :meth:`Molecule.from_rdmol`."""
    mol = MOL.copy()
    mol_rot = molkit.to_rdmol(Molecule(join(PATH, 'Methanol_rotate.xyz')))
    mol.from_rdmol(mol_rot)

    conf = mol_rot.GetConformer()
    for at, at_ref in zip(mol, mol_rot.GetAtoms()):
        pos = conf.GetAtomPosition(at_ref.GetIdx())
        coords = pos.x, pos.y, pos.z
        symbol = at_ref.GetSymbol()
        assertion.eq(at.coords, coords)
        assertion.eq(at.symbol, symbol)


def test_get_index() -> None:
    """Test :meth:`Molecule.get_index`."""
    for j, at in enumerate(MOL, 1):
        i = MOL.get_index(at)
        assertion.eq(i, j)

    ref = [(1, 3), (2, 3), (1, 6), (1, 4), (1, 5)]
    for bond, j in zip(MOL.bonds, ref):
        i = MOL.get_index(bond)
        assertion.eq(i, j)


def test_merge_mol() -> None:
    """Test :meth:`Molecule.merge_mol`."""
    mol = MOL.copy()
    mol_list = [mol.copy() for _ in range(10)]
    atom_list = list(chain.from_iterable(mol_list)) + mol.atoms
    bond_list = list(chain.from_iterable(m.bonds for m in mol_list)) + mol.bonds
    mol.merge_mol(mol_list)

    assertion.eq(len(mol.atoms), len(atom_list))
    assertion.eq(len(mol.bonds), len(bond_list))
    for at in mol.atoms:
        assertion.contains(atom_list, at)
    for bond in mol.bonds:
        assertion.contains(bond_list, bond)


def test_separate_mod() -> None:
    """Test :meth:`Molecule.separate_mod`."""
    mol = MOL.copy()
    mol_list = mol.separate_mod()

    for m in mol_list:
        for at in m.atoms:
            assertion.contains(mol, at)
            assertion.contains(mol.atoms, at)
        for bond in m.bonds:
            assertion.contains(mol.bonds, bond)


def test_to_atnum() -> None:
    """Test :func:`CAT.mol_utils.to_atnum`."""
    for j, (symbol, *_) in enumerate(PeriodicTable.data):
        i = to_atnum(symbol)
        assertion.eq(i, j)
        assertion.eq(to_atnum(j), j)
    assertion.exception(TypeError, to_atnum, {})
    assertion.exception(TypeError, to_atnum, [])
    assertion.exception(TypeError, to_atnum, ())
    assertion.exception(PTError, to_atnum, 'bob')
    assertion.exception(PTError, to_atnum, 'bill')


def test_to_symbol() -> None:
    """Test :func:`CAT.mol_utils.to_symbol`."""
    for j, (symbol, *_) in enumerate(PeriodicTable.data):
        i = to_symbol(j)
        assertion.eq(i, symbol)
        assertion.eq(to_symbol(symbol), symbol)
    assertion.exception(TypeError, to_symbol, {})
    assertion.exception(TypeError, to_symbol, [])
    assertion.exception(TypeError, to_symbol, ())
    assertion.exception(PTError, to_symbol, 999)
    assertion.exception(PTError, to_symbol, -999)


def test_adf_connectivity() -> None:
    """Test :func:`CAT.mol_utils.adf_connectivity`."""
    ref = ['1 3 1.0', '2 3 1.0', '1 6 1.0', '1 4 1.0', '1 5 1.0']
    connectivity_list = adf_connectivity(MOL)
    assertion.eq(connectivity_list, ref)


def test_fix_carboxyl() -> None:
    """Test :func:`CAT.mol_utils.fix_carboxyl`."""
    mol = Molecule(join(PATH, 'Acetate.xyz'))  # Acetate; BP86/QZ4P
    mol.guess_bonds()
    mol[1, 3].order = 2
    mol[1, 4].order = 1
    mol[4].move_to([0.200000, -1.129666, 0.605745])
    mol[4].properties.charge = -1

    fix_carboxyl(mol)
    C, O1, O2 = mol[1], mol[3], mol[4]
    angle = C.angle(O1, O2, result_unit='degree')
    assertion.eq(round(angle), 120)


def test_fix_h() -> None:
    """Test  :func:`CAT.mol_utils.fix_h`."""
    mol = Molecule(join(PATH, 'Ethylene.xyz'))  # Ethylene; BP86/QZ4P
    mol.guess_bonds()
    mol[3].move_to([0.0, 0.3, -0.4])

    fix_h(mol)
    H, C1, C2 = mol[3], mol[1], mol[2]
    angle = C1.angle(H, C2, result_unit='degree')
    assertion.eq(round(angle), 120)
