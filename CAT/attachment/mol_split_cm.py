"""A context manager for temporary splitting (and reassembling) molecules.

Index
-----
.. currentmodule:: CAT.attachment.mol_split_cm
.. autosummary::
    SplitMol

API
---
.. autoclass:: SplitMol
    :members:
    :private-members:
    :special-members:

"""

import copy
import reprlib
from typing import Iterable, Union, Dict, Tuple, NoReturn, Any, Type
from contextlib import AbstractContextManager

from scm.plams import Molecule, Atom, PT, Bond, MoleculeError, PTError, rotation_matrix

from ..mol_utils import separate_mod  # noqa: F401

__all__ = ['SplitMol']


class SplitMol(AbstractContextManager):
    """A context manager for temporary splitting a single molecule into multiple components.

    The context manager splits the provided molecule into multiple components,
    capping all broken bonds in the process.
    The exact amount of fragments depends on the number of specified bonds.

    These moleculair fragments are returned upon opening the context manager and merged back
    into the initial molecule once the context manager is closed. While opened, the initial molecule
    is cleared of all atoms and bonds, while the same hapens to the moleculair fragments upon
    closing.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule, Bond, from_smiles

        >>> mol: Molecule = from_smiles('CC')  # Ethane
        >>> bond: Bond = mol[1, 2]

        # A backup of all bonds and atoms
        >>> bonds_backup = mol.bonds.copy()
        >>> atoms_backup = mol.atoms.copy()

        # The context manager is opened; the bond is removed and the molecule is fragmented
        >>> with SplitMol(mol, bond) as fragment_tuple:
        ...     for fragment in fragment_tuple:
        ...         fancy_operation(fragment)  # doctest: +SKIP
        ...
        ...     print(
        ...         mol.bonds == bonds_backup,
        ...         mol.atoms == atoms_backup,
        ...         bond in mol.bonds
        ...     )
        False False False

        # The context manager is closed; all atoms and bonds have been restored
        >>> print(
        ...     mol.bonds == bonds_backup,
        ...     mol.atoms == atoms_backup,
        ...     bond in mol.bonds
        ... )
        True True True

    Parameters
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.
        See :attr:`SplitMol.mol`.

    bond_list : |plams.Bond| or |Iterable| [|plams.Bond|]
        An iterable consisting of PLAMS bonds.
        All bonds must be part of **mol**.
        See :attr:`SplitMol.bonds`.

    cap_type : :class:`str`, :class:`int` or |plams.Atom|
        An atomic number or symbol of the atom type used for capping the to-be split molecule.
        See :attr:`SplitMol.cap_type`.

    Attributes
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    bonds : :class:`set` [|plams.Bond|]
        A set of PLAMS bonds.

    cap_type : :class:`str`
        An atomic symbol of the atom type used for capping the to-be split molecule.

    _at_pairs : :class:`list` [:class:`dict` [|plams.Atom|, |plams.Atom|]], optional
        A list of dictionaries.
        Each dictionary contains two atoms as keys (see :attr:`SplitMol.bond_list`) and
        their respective capping atom as values.
        Used for reassembling :attr:`SplitMol.mol` once the context manager is closed.
        Set internally by :meth:`SplitMol.__enter__`.

    _vars_backup : :class:`dict` [:class:`str`, |Any|], optional
        A backup of all instance variables of :attr:`SplitMol.mol`.
        Set internally by :meth:`SplitMol.__enter__`.

    _tmp_mol_list : :class:`tuple` [|plams.Molecule|], optional
        A list of PLAMS molecules obtained by splitting :attr:`SplitMol.mol`.
        Set internally by :meth:`SplitMol.__enter__`.

    Raises
    ------
    |MoleculeError|
        Raised when one attempts to access or manipulate the instance variables of
        :attr:`SplitMol.mol` when the context manager is opened.

    """

    def __init__(self, mol: Molecule, bond_list: Union[Bond, Iterable[Bond]],
                 cap_type: Union[str, int, Atom] = 'H') -> None:
        """Initialize a :class:`SplitMol` instance."""
        self.mol = mol
        self.bonds = bond_list
        self.cap_type = cap_type

        # Private variables which are modified in __enter__ and __exit__
        self._at_pairs = None
        self._vars_backup = None
        self._tmp_mol_list = None

    """####################################### Properties #######################################"""

    @property
    def bonds(self) -> Dict[Bond, int]:
        """Getter: Return :attr:`SplitMol.bonds`. Setter: Assign the provided value as a dictionary with bonds as keys and their matching index as value."""  # noqa
        return self._bonds

    @bonds.setter
    def bonds(self, value: Union[Bond, Iterable[Bond]]) -> None:
        self._bonds = {value} if isinstance(value, Bond) else set(value)

    @property
    def cap_type(self) -> str:
        """Getter: Return :attr:`SplitMol.cap_type`. Setter: Assign the provided value after parsing the atom type (:class:`str`)."""  # noqa
        return self._cap_type

    @cap_type.setter
    def cap_type(self, value: Union[str, int, Atom]) -> None:
        if isinstance(value, int):
            self._cap_type = PT.get_symbol(value)

        elif isinstance(value, Atom):
            self._cap_type = value.symbol

        elif isinstance(value, str):
            if value.capitalize() not in PT.symtonum:
                raise PTError(f"Trying to convert an invalid atomic symbol: {reprlib.repr(value)}")
            self._cap_type = value

        else:
            raise TypeError(f"Invalid atom type: {repr(type(value))}")

    """Methods for manipulating the molecule."""

    def cap_fragments(self, bond: Bond) -> Dict[Atom, Atom]:
        """Delete a bond from :attr:`SplitMol.mol` and cap the resulting fragments with :attr:`SplitMol.cap_type`.

        Parameters
        ----------
        bond : |plams.Bond|
            A PLAMS bond.

        Returns
        -------
        :class:`dict` [|plams.Atom|, |plams.Atom|]
            A dictionary with the old atoms in **bond** as keys and their new capping atoms
            as values.

        """  # noqa
        mol = self.mol
        symbol = self.cap_type

        # Construct the capping atoms
        atom1, atom2 = bond.atom1, bond.atom2
        atom1_cap = Atom(symbol=symbol, coords=atom2.coords)
        atom2_cap = Atom(symbol=symbol, coords=atom1.coords)
        mol.add_atom(atom1_cap, adjacent=[atom1])
        mol.add_atom(atom2_cap, adjacent=[atom2])

        # Resize the capping atom bond
        length1 = atom1.radius + atom1_cap.radius
        length2 = atom2.radius + atom2_cap.radius
        mol.bonds[-2].resize(atom1_cap, length1)
        mol.bonds[-1].resize(atom2_cap, length2)

        # Delete the old bond and return a dictionary containg marking all new bonds
        return {atom1: atom1_cap, atom2: atom2_cap}

    def reassemble(self) -> None:
        """Reassemble :attr:`SplitMol.mol` from its consitituent components in :attr:`SplitMol._tmp_mol_list`.

        All capping atoms (and their respective bonds) specified in :attr:`SplitMol._at_pairs` are
        removed and the previously broken bonds (stored in :attr:`SplitMol.bonds`) are restored.

        """  # noqa
        mol = self.mol
        mark = self._at_pairs
        bonds = self.bonds

        for atom_dict, bond in zip(mark, bonds):
            # Extract atoms
            iterator = iter(atom_dict.items())
            atom1, atom1_cap = next(iterator)
            atom2, atom2_cap = next(iterator)

            # Allign the molecules by rotation
            vec1 = atom2.vector_to(atom2_cap)
            vec2 = atom1_cap.vector_to(atom1)
            rotmat = rotation_matrix(vec1, vec2)
            atom2.mol.rotate(rotmat)

            # Allign the molecules by translation
            vec_trans = atom2.vector_to(atom1_cap)
            atom2.mol.translate(vec_trans)

            # Replace the capping atom bonds with the previously broken bond
            atom1.bonds[-1] = bond
            atom2.bonds[-1] = bond

        # Ensure all atoms and bonds belong to mol
        for at in mol.atoms:
            at.mol = mol
        for bond in mol.bonds:
            bond.mol = mol

        # Resize bonds
        for atom_dict, bond in zip(mark, bonds):
            atom1, atom2 = atom_dict
            length = atom1.radius + atom2.radius
            bond.resize(atom1, length)

    def lock_mol(self) -> None:
        """Lock :attr:`SplitMol.mol`, preventing any access to the instance."""
        mol = self.mol
        mol.atoms = mol.bonds = FrozenList()

    def unlock_mol(self) -> None:
        """Unlock :attr:`SplitMol.mol`, restoring access to the instance."""
        self.mol.__dict__ = self._vars_backup

    @staticmethod
    def reset_vars(obj: Any) -> None:
        """Replace all instance variables of **obj** with empty instances of their respective class.

        A value will be substituted for ``None`` if a :exc:`TypeError` is encountered during
        instance creation.

        Parameters
        ----------
        obj : |Any|
            A class instance with the :attr:`__dict__` attribute.

        """
        vars_dct = vars(obj)
        for k, v in vars_dct.items():
            cls = type(v)
            try:  # More foolproof than .__init__(), especially for empty class instances
                vars_dct[k] = cls.__new__(cls)
            except TypeError:
                # For some classes even the .__new__ method fails (e.g. range or np.ndarray objects)
                vars_dct[k] = None

    """############################# Context manager magic methods #############################"""

    def __enter__(self) -> Tuple[Molecule]:
        """Enter the :class:`SplitMol` context manager; return a list of molecules."""
        # Create a backup of mols' instance variables
        mol = self.mol
        self._vars_backup = vars(mol)
        mol.__dict__ = {k: copy.copy(v) for k, v in vars(mol).items()}

        # Add capping atoms along the to-be split bonds
        self._at_pairs = [self.cap_fragments(bond) for bond in self.bonds]

        # Actually delete the to-be split bonds and split the molecule
        for bond in self.bonds:
            mol.delete_bond(bond)
        self._tmp_mol_list = mol.separate_mod()

        # Lock al instance variables of mol
        # Accessing mol will now raise a MoleculeError for as long as the context manager is open
        self.lock_mol()
        return self._tmp_mol_list

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the :class:`SplitMol` context manager."""
        # Restore all instance variables of mol
        self.unlock_mol()

        # Delete all capping atoms recreate the previously broken bond(s)
        self.reassemble()

        # Delete all instance variables of the temporary molecules
        # All variables are replaced with empty instance of their respective class
        for mol_tmp in self._tmp_mol_list:
            self.reset_vars(mol_tmp)


class FrozenList(list):
    """A list subclass that will raise a |MoleculeError| upon modification of an instance."""

    def __init__(self, *args,
                 exc_type: Type[Exception] = MoleculeError,
                 exc_msg: str = ("'Molecule' objects are inaccessible while opened in "
                                 "the 'SplitMol' context manager"),
                 **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.exc_type = exc_type
        self.exc_msg = exc_msg

    def raise_exception(self, *args, **kwargs) -> NoReturn: raise self.exc_type(self.exc_msg)
    __delitem__ = __iadd__ = __setitem__ = raise_exception
    append = clear = extend = insert = pop = remove = reverse = sort = raise_exception
