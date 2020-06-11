"""A context manager for temporary removing a set of atoms from a molecule.

Index
-----
.. currentmodule:: CAT.attachment.remove_atoms_cm
.. autosummary::
    RemoveAtoms

API
---
.. autoclass:: RemoveAtoms
    :members:

"""

import sys
from typing import Iterable, Union, Sequence, cast
from contextlib import AbstractContextManager

from scm.plams import Molecule, Atom, Bond, MoleculeError

if sys.version_info >= (3, 8):
    # dict.__reversed__() was added in Python 3.8
    from builtins import dict as OrderedDict  # noqa: N812
else:
    from collections import OrderedDict

__all__ = ['RemoveAtoms']


class RemoveAtoms(AbstractContextManager):
    """A context manager for temporary removing a set of atoms from a molecule.

    The *relative* ordering of the to-be removed atoms (and matching bonds),
    as specified in **atoms**, is preserved during the removal and reattachment process.
    Note that reattaching will (re-)append the removed atoms/bonds,
    a process which is thus likelly to affect the *absolute* ordering of
    atoms/bonds within the entire molecule.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule, Atom, from_smiles

        >>> mol: Molecule = from_smiles('CO')
        >>> atom1: Atom = mol[1]
        >>> atom2: Atom = mol[2]

        >>> atom_set = {atom1, atom2}
        >>> with RemoveAtoms(mol, atom_set):
        ...     print(atom1 in mol, atom2 in mol)
        False False

        >>> print(atom1 in mol, atom2 in mol)
        True True

    Parameters
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.
        See :attr:`RemoveAtoms.mol`.

    atoms : |plams.Atom| or |Iterable| [|plams.Atom|]
        A PLAMS atom or an iterable consisting of unique PLAMS atoms.
        All supplied atoms should belong to **mol**.
        See :attr:`RemoveAtoms.atoms`.

    Attributes
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    atoms : |Sequence| [|plams.Atom|]
        A sequence of PLAMS atoms belonging to :attr:`RemoveAtoms.mol`.
        Setting a value will convert it into a sequence of atoms.

    _bonds : |OrderedDict| [|plams.Bond|, ``None``]
        A ordered dictionary of PLAMS bonds connected to one or more atoms in
        :attr:`RemoveAtoms.atoms`.
        All values are ``None``, the dictionary serving as an improvised ``OrderedSet``.
        Set to ``None`` until :meth:`RemoveAtoms.__enter__` is called.

    """

    @property
    def atoms(self) -> Sequence[Atom]:
        """Get or set :attr:`RemoveAtoms.atoms`. Setting will convert the supplied value into a sequence of atoms."""  # noqa
        return self._atoms

    @atoms.setter
    def atoms(self, value: Union[Atom, Iterable[Atom]]) -> None:
        if isinstance(value, Atom):  # It's an atom
            self._atoms = (value,)
        elif not hasattr(value, '__reversed__'):  # It's an Iterator or Collection
            self._atoms = tuple(value)
        else:
            self._atoms = value  # It's a Sequence (probably)

    def __init__(self, mol: Molecule, atoms: Union[Atom, Iterable[Atom]]) -> None:
        """Initialize a :class:`RemoveAtoms` instance."""
        self.mol = mol
        self.atoms = cast(Sequence[Atom], atoms)
        self._bonds: 'OrderedDict[Bond, None]' = OrderedDict()  # An improvised "OrderedSet"

    def __enter__(self) -> None:
        """Enter the context manager; delete all atoms in :class:`RemoveAtoms.atoms`."""
        mol = self.mol
        bonds_set = self._bonds
        for atom in self.atoms:
            for bond in atom.bonds:
                bonds_set[bond] = None
            mol.delete_atom(atom)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; reassign all atoms in :class:`RemoveAtoms.atoms`."""
        mol = self.mol
        for atom in reversed(self.atoms):
            mol.add_atom(atom)
        for bond in reversed(self._bonds):
            try:
                mol.add_bond(bond)
            except MoleculeError:
                pass  # One of the bonded atoms has been manually deleted by the user: skip it
        self._bonds = OrderedDict()
