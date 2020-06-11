"""Utility functions related to the :mode:`distribution<CAT.attachment.distribution>` module.

Index
-----
.. currentmodule:: CAT.attachment.distribution_utils
.. autosummary::
    test_distribute

API
---
.. autofunction:: test_distribute

"""

from typing import Any, Union, Optional, Tuple, Iterable
from collections import abc

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule, Atom, rotation_matrix

from .distribution import distribute_idx


def test_distribute(mol: Union[Molecule, str], symbol: str,
                    f_range: Union[float, Iterable[float]],
                    rotate: Optional[Tuple[float, float, float]] = (0.1, -0.1, 0.9),
                    **kwargs: Any) -> Molecule:
    r"""Test function for :func:`CAT.attachment.distribution.distribute_idx`.

    Examples
    --------
    .. code:: python

        >>> import numpy as np
        >>> from scm.plams import Molecule

        >>> mol_input = Molecule()
        >>> xyz_output: str = ...
        >>> at_symbol = 'Cl'
        >>> f_range: np.ndarray = 2**-np.arange(8.0)

        >>> mol_out: Molecule = test_distribute(mol_input, at_symbol, f_range)  # doctest: +SKIP
        >>> mol_out.write(xyz_output)  # doctest: +SKIP

        >>> print(len(mol_input) == len(p_range) * len(mol_out))  # doctest: +SKIP
        True

    Parameters
    ----------
    mol : :class:`Molecule` or :class:`str`
        A molecule or path+filename containing a molecule.

    symbol : :class:`str`
        The atomic symbol of the anchor atom.

    f_range : :class:`float` or :class:`Iterable<collections.abc.Iterable>` [:class:`float`]
        A float or iterable of floats subject to the following constraint: :math:`0 < f \le 1`.

    rotate : :class:`Sequence<collection.abc.Sequence>` [:class:`float`], shape :math:`(3,)`, optional
        A sequence of three floats representing a molecular orientation.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for
        :func:`distribute_idx()<CAT.attachment.distribution.distribute_idx>`:
        ``follow_edge``, ``mode`` and ``start``.

    Returns
    -------
    :class:`Molecule`
        A Molecule instance containing one molecule for every item in **p_range**

    """  # noqa
    if not isinstance(mol, Molecule):
        mol = Molecule(mol)
    if not isinstance(f_range, abc.Iterable):
        f_range = (f_range,)

    ret = Molecule()
    trans = cdist(mol, mol).max() * 1.1
    for i, f in enumerate(f_range):
        mol_tmp = _test_distribute(mol, symbol, f=f, **kwargs)
        if rotate is not None:
            mol_tmp.rotate(rotation_matrix([0, 0, 1], rotate))
        mol_tmp.translate([i*trans, 0, 0])
        ret += mol_tmp
    return ret


def _test_distribute(mol: Molecule, symbol: str, **kwargs) -> Molecule:
    """Helper function for :func:`test_distribute`."""
    if not isinstance(mol, Molecule):
        mol = Molecule(mol)

    _idx_in = [i for i, at in enumerate(mol) if at.symbol == symbol]
    idx_in = np.fromiter(_idx_in, count=len(_idx_in), dtype=int)
    idx_out = distribute_idx(mol, idx_in, **kwargs)

    a = symbol
    b = 'I' if a != 'I' else 'Br'
    mol2 = Molecule()
    for i, at in enumerate(mol):
        if at.symbol != symbol:
            continue
        symbol_new = a if i not in idx_out else b
        mol2.add_atom(Atom(symbol=symbol_new, coords=at.coords, mol=mol2))
    return mol2
