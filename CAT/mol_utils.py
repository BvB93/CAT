"""A module with misc functions related to manipulating molecules and their geometry.

Index
-----
.. currentmodule:: CAT.mol_utils
.. autosummary::
    from_mol_other
    from_rdmol
    get_index
    separate_mod
    to_atnum
    to_symbol
    adf_connectivity
    fix_carboxyl
    fix_h

API
---
.. automethod:: from_mol_other
.. automethod:: from_rdmol
.. automethod:: get_index
.. automethod:: separate_mod
.. autofunction:: to_atnum
.. autofunction:: to_symbol
.. autofunction:: adf_connectivity
.. autofunction:: fix_carboxyl
.. autofunction:: fix_h

"""

from typing import Optional, Iterable, Union, Tuple, List

import numpy as np

from scm.plams import Molecule, Atom, Bond, MoleculeError, add_to_class, PTError
from scm.plams.tools.periodic_table import PeriodicTable
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

__all__ = ['adf_connectivity', 'fix_h', 'fix_carboxyl']

try:
    SANITIZE: int = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
except TypeError:  # This prevents Sphinx from raising a TypeError when rdkit is mocked
    SANITIZE = 0


@add_to_class(Molecule)
def from_mol_other(self, mol: Molecule,
                   atom_subset: Optional[Iterable[Atom]] = None) -> None:
    """Update the Cartesian coordinates of this instance with those from another PLAMS molecule.

    Alternatively, update only a subset of atoms.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    atom_subset : |list|_ [|plams.Atom|_]
        Optional: A subset of atoms in **self**.

    """
    at_subset = atom_subset or self.atoms
    for at1, at2 in zip(at_subset, mol):
        at1.coords = at2.coords


@add_to_class(Molecule)
def from_rdmol(self, rdmol: Chem.Mol,
               atom_subset: Optional[Iterable[Atom]] = None) -> None:
    """Update the atomic coordinates of this instance with coordinates from an RDKit molecule.

    Alternatively, update only a subset of atoms.

    Parameters
    ----------
    rdmol : |rdkit.Chem.Mol|_
        An RDKit molecule.

    atom_subset : |list|_ [|plams.Atom|_]
        Optional: A subset of atoms in **self**.

    """
    at_subset = atom_subset or self.atoms
    conf = rdmol.GetConformer()
    for at1, at2 in zip(at_subset, rdmol.GetAtoms()):
        pos = conf.GetAtomPosition(at2.GetIdx())
        at1.coords = (pos.x, pos.y, pos.z)


@add_to_class(Molecule)
def get_index(self, value: Union[Atom, Bond]) -> Union[int, Tuple[int, int]]:
    """Return the first index of **value** within this instance.

    **value** expects an instance of either :class:`Atom` or :class:`Bond`.

    Note
    ----
    Following the convention addopted by PLAMS, the returned index/indices are 1-based rather
    than 0-based.

    Parameters
    ----------
    value : |plams.Atom|_ or |plams.Bond|_
        A PLAMS atom or bonds.

    Returns
    -------
    |int|_ or |tuple|_ [|int|_]
        An atomic index or (**value**: |plams.Atom|_) or
        a tuple of two atomic indices (**item**: |plams.Bond|_).

    Raises
    ------
    TypeError
        Raised if **value** is an instance of neither :class:`Atom` nor :class:`Bond`.

    MoleculeError
        Raised if the passed :class:`Atom` or :class:`Bond` is not in this instance.

    """
    if isinstance(value, Atom):
        if value not in self.atoms:
            raise MoleculeError("Passed atom, {repr(value)}, is not in this instance")
        return 1 + self.atoms.index(value)
    elif isinstance(value, Bond):
        if value not in self.bonds:
            raise MoleculeError(f"Passed bond, {repr(value)}, is not in this instance")
        at1, at2 = value
        return 1 + self.atoms.index(at1), 1 + self.atoms.index(at2)

    err = "item excepts an instance of 'Atom' or 'Bond'; observed type: '{}'"
    raise TypeError(err.format(value.__class__.__name__))


@add_to_class(Molecule)
def separate_mod(self) -> Tuple[Molecule]:
    """Modified PLAMS function: creates new molecules out of this instance rather than a copy of this instance. Atoms, bonds and properties are *not* copied.

    Separate the molecule into connected components.
    Returns is a list of new Molecule instrances (all atoms and bonds are disjoint with
    the original molecule).
    Each element of this list is identical to one connected component of the base molecule.
    A connected component is a subset of atoms such that there exists a path
    (along one or more bonds) between any two atoms.

    Returns
    -------
    |tuple|_ [|plams.Molecule|_]
        A list of molecules with atoms and bonds from **self**.

    """  # noqa: E501
    frags = ()
    for at in self:
        at._visited = False

    def dfs(v, mol):
        v._visited = True
        v.mol = mol
        for e in v.bonds:
            e.mol = mol
            u = e.other_end(v)
            if not u._visited:
                dfs(u, mol)

    for src in self.atoms:
        if not src._visited:
            m = Molecule()
            dfs(src, m)
            frags += (m,)
            m.properties = self.properties.copy()

    for at in self.atoms:
        del at._visited
        at.mol.atoms.append(at)
    for b in self.bonds:
        b.mol.bonds.append(b)

    return frags


@add_to_class(Molecule)
def round_coords(self, decimals: int = 3) -> None:
    """Round the Cartesian coordinates of this instance to a given precision in decimal digits.

    Performs an inplace update of all atoms in this instance.

    Parameters
    ----------
    decimals : int
        The desired precision in decimal digits.

    """
    xyz = self.as_array()
    xyz[:] = np.round(xyz, decimals=decimals)
    self.from_array(xyz)


def to_atnum(value: Union[str, int, Atom]) -> int:
    """Turn an atomic symbol into an atomic number.

    Parameters
    ----------
    value : |int|_ or |str|_
        An Atom, atomic symbol or number.

    Returns
    -------
    |int|_
        An atomic number.

    Raises
    ------
    TypeError
        Raised if **item** is an instance of neither :class:`str` nor :class:`int`.

    """
    if isinstance(value, str):
        try:
            return PeriodicTable.get_atomic_number(value)
        except PTError as ex:
            ex.args = (ex.args[0] + f': {repr(value)}',)
            raise ex
    elif isinstance(value, int):
        if value not in PeriodicTable.symtonum.values():
            raise PTError(f'trying to convert incorrect atomic number: {repr(value)}')
        return value

    try:
        return value.atnum
    except AttributeError as ex:
        raise TypeError("'value' expected an Atom or atomic symbol/number; "
                        f"observed type: '{value.__class__.__name__}'") from ex


def to_symbol(value: Union[str, int, Atom]) -> str:
    """Turn an atomic number into an atomic symbol.

    Parameters
    ----------
    value : |int|_ or |str|_
        An Atom, atomic symbol or number.

    Returns
    -------
    |int|_
        An atomic symbol.

    Raises
    ------
    TypeError
        Raised if **item** is an instance of neither :class:`str` nor :class:`int`.

    """
    if isinstance(value, int):
        try:
            return PeriodicTable.get_symbol(value)
        except PTError as ex:
            ex.args = (ex.args[0] + f': {repr(value)}',)
            raise ex
    elif isinstance(value, str):
        if value not in PeriodicTable.symtonum.keys():
            raise PTError(f'trying to convert incorrect atomic symbol: {repr(value)}')
        return value

    try:
        return value.symbol
    except AttributeError as ex:
        raise TypeError("'value' expected an Atom or atomic symbol/number; "
                        f"observed type: '{value.__class__.__name__}'") from ex


def adf_connectivity(mol: Molecule) -> List[str]:
    """Create an AMS-compatible connectivity list.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule with :math:`n` bonds.

    Returns
    -------
    :math:`n` |list|_ [|str|_]
        An ADF-compatible connectivity list of :math:`n` bonds.

    """
    mol.set_atoms_id()

    # Create list of indices of all aromatic bonds
    try:
        rdmol = molkit.to_rdmol(mol)
    except Exception as ex:
        if type(ex) is ValueError or ex.__class__.__name__ == 'ArgumentError':
            # Plan B: ignore aromatic bonds
            bonds = [f'{bond.atom1.id} {bond.atom2.id} {bond.order:.1f}' for bond in mol.bonds]
            mol.unset_atoms_id()
            return bonds
        raise ex

    aromatic = [bond.GetIsAromatic() for bond in rdmol.GetBonds()]

    # Create a list of bond orders; aromatic bonds get a bond order of 1.5
    bond_orders = [(1.5 if ar else bond.order) for ar, bond in zip(aromatic, mol.bonds)]
    bonds = [f'{bond.atom1.id} {bond.atom2.id} {order:.1f}' for
             bond, order in zip(mol.bonds, bond_orders)]
    mol.unset_atoms_id()

    return bonds


def _smiles_to_rdmol(smiles: str) -> Chem.Mol:
    """Convert a SMILES string into an rdkit Mol; supports explicit hydrogens."""
    # RDKit tends to remove explicit hydrogens if SANITIZE_ADJUSTHS is enabled
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.rdmolops.SanitizeMol(mol, sanitizeOps=SANITIZE)
    except Exception as ex:
        raise ex.__class__(f'Failed to parse the following SMILES string: {repr(smiles)}\n\n{ex}')
    return mol


#: A carboxylate
_CARBOXYLATE: Chem.Mol = _smiles_to_rdmol('[O-]C(C)=O')


def fix_carboxyl(mol: Molecule) -> None:
    """Resets carboxylate OCO angles if it is smaller than :math:`60` degrees.

    Performs an inplace update of **plams_mol**.

    Parameters
    ----------
    plams_mol : |plams.Molecule|_
        A PLAMS molecule.

    """
    rdmol = molkit.to_rdmol(mol)
    conf = rdmol.GetConformer()
    matches = rdmol.GetSubstructMatches(_CARBOXYLATE)

    if matches:
        get_angle = rdMolTransforms.GetAngleDeg
        set_angle = rdMolTransforms.SetAngleDeg
        for idx in matches:
            if get_angle(conf, idx[3], idx[1], idx[0]) < 60:
                set_angle(conf, idx[2], idx[1], idx[3], 180.0)
                set_angle(conf, idx[0], idx[1], idx[3], 120.0)
        mol.from_rdmol(rdmol)


def fix_h(mol: Molecule) -> None:
    """If a C=C-H angle is smaller than :math:`20` degrees, set it back to :math:`120` degrees.

    Performs an inplace update of **plams_mol**.

    Parameters
    ----------
    plams_mol : |plams.Molecule|_
        A PLAMS molecule.

    """
    h_list = [atom for atom in mol if atom.atnum == 1 and 2.0 in
              [bond.order for bond in mol.neighbors(atom)[0].bonds]]

    rdmol = molkit.to_rdmol(mol)
    conf = rdmol.GetConformer()
    get_idx = mol.atoms.index
    set_angle = rdMolTransforms.SetAngleDeg
    get_angle = rdMolTransforms.GetAngleDeg

    update = False
    for atom in h_list:
        at1 = atom  # Central atom
        at2 = mol.neighbors(at1)[0]  # Neighbours
        at3 = [atom for atom in mol.neighbors(at2) if atom != at1]  # Neighbours of neighbours

        # Create 2 sets of 3 atomic indices for defining angles: at1-at2=at3
        idx_tup1 = get_idx(at3[0]), get_idx(at2), get_idx(at1)
        idx_tup2 = get_idx(at3[1]), get_idx(at2), get_idx(at1)

        if get_angle(conf, *idx_tup1) <= 20.0:
            set_angle(conf, *idx_tup1, 120.0)
            update = True
        elif get_angle(conf, *idx_tup2) <= 20.0:
            set_angle(conf, *idx_tup2, 120.0)
            update = True

    if update:
        mol.from_rdmol(rdmol)
