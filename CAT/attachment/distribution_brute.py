"""Functions for creating distributions of atomic indices using brute-force approaches.

Index
-----
.. currentmodule:: CAT.attachment.distribution_brute
.. autosummary::
    brute_uniform_idx

API
---
.. autofunction:: brute_uniform_idx

"""

import reprlib
from typing import Union, Callable

import numpy as np

from scm.plams import Molecule
from nanoutils import array_combinations

from .distribution import OPERATION_MAPPING

__all__ = ['brute_uniform_idx']


def brute_uniform_idx(mol: Union[Molecule, np.ndarray],
                      idx: Union[int, np.ndarray], n: int = 2,
                      operation: str = 'min',
                      weight: Callable[[np.ndarray], np.ndarray] = lambda x: np.exp(-x)
                      ) -> np.ndarray:
    r"""Brute force approach to creating uniform or clustered distributions.

    Explores, and evaluates, all valid combinations of size :math:`n` constructed
    from the :math:`k` atoms in **neighbor** closest to each atom in **center**.

    The combination where the :math:`n` atoms are closest (:code:`operation = 'max'`) or
    furthest removed from each other (:code:`operation = 'min'`) is returned.

    Parameters
    ----------
    mol : array-like [:class:`float`], shape :math:`(m,3)`
        An array-like object with Cartesian coordinate representing a collection of central atoms.

    idx : array-like [:class:`int`], shape :math:`(l,p)`
        An array-like object with indices in **mol**.
        Combinations will be explored and evaluated along axis ``-1`` of the passed array.

    n : :class:`int`
        The number of to-be returned opposing atoms.
        Should be larger than or equal to 1.

    operation : :class:`str`
        Whether to evaluate the weighted distance using :func:`argmin()<numpy.nanargmin>` or
        :func:`argmax()<numpy.nanargmax>`.
        Accepted values are ``"min"`` and ``"max"``.

    weight : :data:`Callable<typing.Callable>`
        A callable for applying weights to the distance; default: :math:`e^{-x}`.
        The callable should take an array as argument and return a new array,
        *e.g.* :func:`numpy.exp`.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`int`], shape :math:`(m, n)`
        An array with indices of opposing atoms.

    See Also
    --------
    :func:`uniform_idx()<CAT.attachment.distribution.uniform_idx>`
        Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    """  # noqa
    try:  # Parse **operation**
        arg_func = OPERATION_MAPPING[operation]
    except KeyError as ex:
        raise ValueError(f"Invalid value for 'operation' ({reprlib.repr(operation)}); "
                         "accepted values: ('min', 'max')") from ex

    # Find the n atoms in mol2 closest to each atom in mol1
    idx = np.array(idx, ndmin=1, copy=False)
    xyz = np.array(mol, ndmin=2, copy=False)
    if not (0 < n <= idx.shape[-1]):
        raise ValueError("'n' should be larger than 0 and smaller than or equal to the last axis of"
                         f" 'idx' ({repr(idx.shape[-1])}); observed value: {repr(n)}")

    # Evaluate all combinations of length n constructed from an iterable of size k
    idx2 = np.swapaxes(array_combinations(idx, r=n), 0, 1)
    xyz2 = xyz[idx2]

    # Construct a weighted distance matrix
    dist = np.triu(np.linalg.norm(xyz2[..., None, :] - xyz2[..., None, :, :], axis=-1))
    dist.shape = len(idx), -1, n**2
    dist_weight = weight(dist).sum(axis=-1)

    # Slice and return
    i = np.arange(len(idx))
    j = arg_func(dist_weight, axis=-1)
    return idx2[i, j]
