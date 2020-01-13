"""
CAT.attachment.distribution
===========================

Functions for creating distributions of atomic indices (*i.e.* core anchor atoms).

Index
-----
.. currentmodule:: CAT.attachment.distribution
.. autosummary::
    distribute_idx
    uniform_idx
    random_idx

API
---
.. autofunction:: distribute_idx
.. autofunction:: uniform_idx
.. autofunction:: random_idx

"""

import reprlib
from typing import Generator, Optional, Iterable, FrozenSet, Any, Union, Callable
from itertools import islice
from collections import abc

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule, Atom

from CAT.attachment.edge_distance import edge_dist

__all__ = ['distribute_idx']

#: A set of allowed values for the **mode** parameter in :func:`get_distribution`.
MODE_SET: FrozenSet[str] = frozenset({'uniform', 'random', 'cluster'})


def distribute_idx(core: Union[Molecule, np.ndarray], idx: Union[int, Iterable[int]], p: float,
                   mode: str = 'uniform', **kwargs: Any) -> np.ndarray:
    r"""Create a new distribution of atomic indices from **idx** of length :code:`p * len(idx)`.

    Parameters
    ----------
    core : :math:`(m, 3)` array-like [:class:`float`]
        A 2D array-like object (such as a :class:`Molecule` instance) consisting
        of Cartesian coordinates.

    idx : :class:`int` or :math:`(i,)` :class:`Iterable<collections.abc.Iterable>` [:class:`int`]
        An integer or iterable of unique integers representing the 0-based indices of
        all anchor atoms in **core**.

    p : :class:`float`
        A float obeying the following condition: :math:`0.0 < p <= 1.0`.
        Represents the fraction of **idx** that will be returned.

    mode : :class:`str`
        How the subset of to-be returned indices will be generated.
        Accepts one of the following values:

        * ``"random"``: A random distribution.
        * ``"uniform"``: A uniform distribution; the distance between each successive atom and
          all previous points is maximized.
        * ``"cluster"``: A clustered distribution; the distance between each successive atom and
          all previous points is minmized.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for the **mode**-specific functions.

    Returns
    -------
    :math:`(p*i,)` :class:`numpy.ndarray` [:class:`int`]
        A 1D array of atomic indices.
        If **idx** has :math:`i` elements,
        then the length of the returned list is equal to :math:`\max(1, p*i)`.

    See Also
    --------
    :func:`uniform_idx`
        Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    :func:`cluster_idx`
        Return the column-indices of **dist** which yield a clustered distribution.

    """
    # Convert **idx** into an array
    try:
        idx_ar = np.array(idx, dtype=int, ndmin=1, copy=False)
    except TypeError:  # A Collection or Iterator
        idx_ar = np.fromiter(idx, dtype=int)

    # Validate the input
    if mode not in MODE_SET:
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(mode)}); "
                         f"accepted values: {reprlib.repr(tuple(MODE_SET))}")
    elif not (0.0 < p <= 1.0):
        raise ValueError("'p' should be larger than 0.0 and smaller than or equal to 1.0; "
                         f"observed value: {reprlib.repr(p)}")
    elif p == 1.0:  # Ensure that **idx** is always returned as copy
        return idx_ar.copy() if idx_ar is idx else idx_ar

    # Create an array of indices
    stop = max(1, int(round(p * len(idx_ar))))
    if mode in ('uniform', 'cluster'):
        xyz = np.array(core, dtype=float, ndmin=2, copy=False)[idx_ar]
        dist = edge_dist(xyz) if kwargs.get('follow_edge', False) else cdist(xyz, xyz)
        operation = 'max' if mode == 'uniform' else 'min'
        generator1 = uniform_idx(dist, operation=operation,
                                 start=kwargs.get('start', None),
                                 cluster_size=kwargs.get('cluster_size', 1))
        generator2 = islice(generator1, stop)
        ret = idx_ar[np.fromiter(generator2, count=stop, dtype=int)]

    elif mode == 'random':
        ret = np.random.permutation(idx_ar)

    # Return a list of `p * len(idx)` atomic indices
    return ret[:stop]


def uniform_idx(dist: np.ndarray, operation: str = 'max', p: float = -2.0,
                cluster_size: int = 1, start: Optional[int] = None) -> Generator[int, None, None]:
    r"""Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    Given the (symmetric) distance matrix :math:`\boldsymbol{D} \in \mathbb{R}^{n,n}` and
    the vector :math:`\boldsymbol{d} \in \mathbb{N}^{\le n}`
    (representing a subset of indices in :math:`D`),
    then the :math:`i`'th element :math:`\boldsymbol{d}_{i}` is
    defined as following:

    .. math::

        \DeclareMathOperator*{\argmax}{\arg\!\max}
        d_{i} = \begin{cases}
            \argmax\limits_{k \in \mathbb{N}} || \boldsymbol{D}_{k,:} ||_{p} &&&
            \text{if} & i=0 \\
            \argmax\limits_{k \in \mathbb{N}} || \boldsymbol{D}[k; d_{0},...,d_{i-1}] ||_{p} &
            \text{with} & k \notin \boldsymbol{d}[0, ..., i-1] &
            \text{if} & i > 0, {i \over m} \in \mathbb{N} \\
            \argmax\limits_{k \in \mathbb{N}}
                || \boldsymbol{D}[k; d_{0},...,d_{i-m}] ||_{p} *
                || \boldsymbol{D}[k; d_{i-m+1},...,d_{i-1}] ||_{p}^{-1} &
            \text{with} & k \notin \boldsymbol{d}[0, ..., i-1] &
            \text{if} & i > 0, {i \over m} \notin \mathbb{Z}
        \end{cases}

    By default :math:`p=-2`.
    Using a negative Minkowski norm is equivalent to, temporarily, projecting the distance matrix
    into recipropal space, thus results in an increased weight of all neighbouring atoms.

    The row in :math:`D` corresponding to :math:`d_{0}`
    can alternatively be specified by **start**.

    The :math:`\text{argmax}` operation can be exchanged for :math:`\text{argmin}` by settings
    **operation** to ``"min"``, thus yielding a clustered- rather than uniform-distribution.

    Parameters
    ----------
    dist : :math:`(m, m)` :class:`numpy.ndarray` [:class:`float`]
        A symmetric 2D NumPy array (:code:`(dist == dist.T).all()`) representing the
        distance matrix :math:`D`.

    operation : :class:`str`
        Whether to minimize or maximize the distance between points.
        Accepted values are ``"min"`` and ``"max"``.

    start : :class:`int`, optional
        The index of the starting row in **dist**.
        If ``None``, start in whichever row contains the global minimum
        (:math:`\DeclareMathOperator*{\argmin}{\arg\!\min} \argmin\limits_{k \in \mathbb{N}} ||\boldsymbol{D}_{k, :}||`) or maximum
        (:math:`\DeclareMathOperator*{\argmax}{\arg\!\max} \argmax\limits_{k \in \mathbb{N}} ||\boldsymbol{D}_{k, :}||`).
        See **operation**.

    p : :class:`float`
        The order of the Minkowski norm; used for determining the optimal values of :math:`d_{i>0}`.
        :math:`p=2` is equivalent to the Euclidian norm:

        .. math::

            || \boldsymbol{x} ||_{p} = \left( \sum_{i=0}^n {x_{i}}^{p} \right)^{1/p}
            \quad \text{with} \quad \boldsymbol{x} \in \mathbb{R}^n

    Yields
    ------
    :class:`int`
        Column-indices specified in :math:`\boldsymbol{d}`.

    """  # noqa
    if operation not in ('min', 'max'):
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(operation)}); "
                         f"accepted values: ('min', 'max')")
    p_inv = 1 / p

    # Truncate and square the distance matrix
    dist_sqr = np.array(dist, dtype=float, copy=True)
    np.fill_diagonal(dist_sqr, np.nan)
    dist_sqr **= p

    # Use either argmin or argmax
    arg_func = np.nanargin if operation == 'min' else np.nanargmax
    start = arg_func(np.nansum(dist_sqr, axis=1)**p_inv) if start is None else start

    # Yield the first index
    dist_1d_sqr = dist_sqr[start].copy()
    dist_1d_sqr[start] = np.nan
    yield start

    # Construct a generator for yielding the remaining indices
    if cluster_size == 1:
        generator = _min_or_max(dist_sqr, dist_1d_sqr, arg_func, p_inv)
    else:
        generator = _min_and_max(dist_sqr, dist_1d_sqr, arg_func, p_inv, cluster_size)

    # Yield remaining indices
    for i in generator:
        yield i


def _min_or_max(dist_sqr: np.ndarray, dist_1d_sqr: np.ndarray,
                arg_func: Callable[[np.ndarray], int], p_inv: float = -0.5
                ) -> Generator[int, None, None]:
    """Helper function for :func:`uniform_idx` if :code:`cluster_size == 1`."""
    for _ in range(len(dist_1d_sqr)-1):
        dist_1d = dist_1d_sqr**p_inv
        i = arg_func(dist_1d)
        dist_1d_sqr[i] = np.nan
        dist_1d_sqr += dist_sqr[i]
        yield i


def _min_and_max(dist_sqr: np.ndarray, dist_1d_sqr: np.ndarray,
                 arg_func: Callable[[np.ndarray], int], p_inv: float = -0.5,
                 cluster_size: int = 1) -> Generator[int, None, None]:
    """Helper function for :func:`uniform_idx` if :code:`cluster_size != 1`."""
    bool_ar = np.zeros(len(dist_1d_sqr)-1, dtype=bool)
    bool_ar[::cluster_size] = True
    j_ar = np.zeros(len(dist_1d_sqr), dtype=float)
    for i in bool_ar:
        if i:
            dist_1d_sqr += j_ar
            j_ar[:] = 0.0
            dist_1d = dist_1d_sqr**p_inv
        else:
            dist_1d = dist_1d_sqr**p_inv
            dist_1d /= j_ar**p_inv
        j = arg_func(dist_1d)
        dist_1d_sqr[j] = np.nan
        j_ar += dist_sqr[j]
        yield j


def cluster_idx(dist: np.ndarray, start: Optional[int] = None) -> np.ndarray:
    r"""Return the column-indices of **dist** which yield a clustered distribution.

    Given the (symmetric) distance matrix :math:`D \in \mathbb{R}^{n, n}` and the starting row
    :math:`\DeclareMathOperator*{\argmin}{\arg\!\min} i = \argmin_{i} ||D_{i, :}||_{2}`,
    return the column-indices of :math:`D_{i, :}` sorted in order of ascending distance.

    .. math::

        \DeclareMathOperator*{\argmin}{\arg\!\min}
        d_{i} = \begin{cases}
            \argmin\limits_{k \in \mathbb{N}} || \boldsymbol{D}_{k,:}|| &&&
            \text{if} & i=0 \\
            \argmin\limits_{k \in \mathbb{N}} D_{k, d_{0}} &
            \text{with} & k \notin \boldsymbol{d}[0, ..., i-1] &
            \text{if} & i \ne 0
        \end{cases}

    Parameters
    ----------
    dist : :math:`(m, m)` array-like [:class:`float`]
        A symmetric 2D NumPy array (:code:`(dist == dist.T).all()`) representing the
        distance matrix :math:`D`.

    start : :class:`int`, optional
        The index of the starting row in **dist**.
        If ``None``, start in row:
        :math:`\DeclareMathOperator*{\argmin}{\arg\!\min} \argmin_{k \in \mathbb{N}} ||D_{k, :}||`.

    Returns
    -------
    :math:`(m,)` :class:`numpy.ndarray` [:class:`int`]
        A 1D array of indices.

    """
    dist = np.asarray(dist, dtype=float)
    start = np.linalg.norm(dist, axis=1).argmin() if start is None else start

    r = dist[start]
    r_arg = r.argsort()
    idx = np.arange(len(r))
    return idx[r_arg]


def _test_distribute(mol: Molecule, symbol: str, **kwargs) -> Molecule:
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


def test_distribute(mol: Union[Molecule, str], symbol: str,
                    p_range: Union[float, Iterable[float]], **kwargs) -> Molecule:
    r"""Test function for :func:`CAT.attachment.distribution.distribute_idx`.

    Examples
    --------
    .. code:: python

        >>> import numpy as np
        >>> from scm.plams import Molecule

        >>> mol_input: Molecule = Molecule(...)
        >>> xyz_output: str = ...
        >>> at_symbol: str = 'Cl'
        >>> p_range: numpy.ndarray = 2**-np.arange(8.0)

        >>> mol_out: Molecule = test_distribute(mol_input, at_symbol, p_range)
        >>> mol_out.write(xyz_output)

        >>> print(len(mol_input) == len(p_range) * len(mol_out))
        True

    Parameters
    ----------
    mol : :class:`Molecule` or :class:`str`
        A molecule or path+filename containing a molecule.

    symbol : :class:`str`
        The atomic symbol of the anchor atom.

    p_range : :class:`float` or :class:`Iterable<collections.abc.Iterable>` :class:`float`
        A float or iterable of floats subject to the following constraint: :math:`0 < p \le 1`.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for :func:`CAT.attachment.distribution.distribute_idx`:
        ``follow_edge``, ``mode`` and ``start``.

    Returns
    -------
    :class:`Molecule`
        A Molecule instance containing one molecule for every item in **p_range**

    """
    if not isinstance(mol, Molecule):
        mol = Molecule(mol)
    if not isinstance(p_range, abc.Iterable):
        p_range = (p_range,)

    ret = Molecule()
    trans = cdist(mol, mol).max() * 1.1
    for i, p in enumerate(p_range):
        mol_tmp = _test_distribute(mol, symbol, p=p, **kwargs)
        mol_tmp.translate([i*trans, 0, 0])
        ret += mol_tmp
    return ret


file = r"/Users/bvanbeek/Downloads/8nm_model_cb_withdummy(1).xyz"
p_range = 2**-np.arange(1.0, 5.0)
mol = test_distribute(file, 'Cl', p_range=p_range, follow_edge=True, mode='uniform', cluster_size=4)
mol.from_array(mol.as_array() / 2)
mol.write(file.replace('xyz', 'output.xyz'))
