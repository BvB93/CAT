"""A context manager for temporary interconverting between PLAMS molecules and NumPy arrays.

Index
-----
.. currentmodule:: CAT.attachment.as_array
.. autosummary::
    AsArray

API
---
.. autoclass:: AsArray
    :members:
    :private-members:
    :special-members:

"""

from typing import Iterable, Union, Sequence
from contextlib import AbstractContextManager
from collections import abc

import numpy as np

from scm.plams import Atom, Molecule

__all__ = ['AsArray']


class AsArray(AbstractContextManager):
    r"""A context manager for temporary interconverting between PLAMS molecules and NumPy arrays.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule

        # Create a H2 example molecule
        >>> h1 = Atom(symbol='H', coords=(0.0, 0.0, 0.0))
        >>> h2 = Atom(symbol='H', coords=(1.0, 0.0, 0.0))
        >>> mol = Molecule()
        >>> mol.add_atom(h1)
        >>> mol.add_atom(h2)

        >>> print(mol)  # doctest: +SKIP
          Atoms:
            1         H      0.000000      0.000000      0.000000
            2         H      1.000000      0.000000      0.000000
        <BLANKLINE>

        # Example: Translate the molecule along the Cartesian Z-axis by 5 Angstroem
        >>> with AsArray(mol) as xyz:
        ...     xyz[:, 2] += 5

        >>> print(mol)  # doctest: +SKIP
          Atoms:
            1         H      0.000000      0.000000      5.000000
            2         H      1.000000      0.000000      5.000000
        <BLANKLINE>

    Parameters
    ----------
    mol : |plams.Molecule| or |Iterable| [|plams.Atom|]
        An iterable consisting of PLAMS atoms.
        See :attr:`AsArray.mol`.

    Attributes
    ----------
    mol : |plams.Molecule| or |Sequence| [|plams.Atom|]
        A PLAMS molecule or a sequence of PLAMS atoms.

    _xyz : :math:`n*3` :class:`numpy.ndarray` [:class:`float`], optional
        A 2D array with the Cartesian coordinates of **mol**.
        Empty by default; this value is set internally by the :meth:`AsArray.__enter__` method.

    """

    @property
    def mol(self) -> Union[Molecule, Sequence[Atom]]:
        """Get or set the embedded molecule."""
        return self._mol

    @mol.setter
    def mol(self, value: Iterable[Atom]) -> None:
        self._mol = value if isinstance(value, (abc.Sequence, Molecule)) else tuple(value)

    def __init__(self, mol: Iterable[Atom]) -> None:
        """Initialize a :class:`AsArray` instance."""
        self.mol = mol
        self._xyz = None

    def __enter__(self) -> np.ndarray:
        """Enter the context manager; return an array of Cartesian coordinates."""
        self._xyz = Molecule.as_array(None, atom_subset=self.mol)
        return self._xyz

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager; update the Cartesian coordinates of :attr:`AsArray.mol`."""
        Molecule.from_array(None, self._xyz, atom_subset=self.mol)
        self._xyz = None
