"""Functions for creating distributions of atomic indices (*i.e.* core anchor atoms).

Index
-----
.. currentmodule:: CAT.attachment.distribution
.. autosummary::
    distribute_idx
    uniform_idx

API
---
.. autofunction:: distribute_idx
.. autofunction:: uniform_idx

"""

import reprlib
import functools
from types import MappingProxyType
from typing import Generator, Optional, Iterable, FrozenSet, Any, Union, Callable, Mapping, TYPE_CHECKING
from itertools import islice, takewhile
from collections import abc

import numpy as np
from scipy.spatial.distance import cdist

from scm.plams import Molecule
from nanoutils import as_nd_array, Literal

from .edge_distance import edge_dist
from ..utils import cycle_accumulate

if TYPE_CHECKING:
    from numpy.typing import ArrayLike
else:
    ArrayLike = 'numpy.typing.ArrayLike'

__all__ = ['distribute_idx']

#: A set of allowed values for the **mode** parameter in :func:`get_distribution`.
MODE_SET: FrozenSet[str] = frozenset({'uniform', 'random', 'cluster'})

#: Map the **operation** parameter in :func:`uniform_idx` to either
#: :func:`numpy.nanargmin` or :func:`numpy.nanargmax`.
OPERATION_MAPPING: Mapping[Union[str, Callable], Callable] = MappingProxyType({
    'min': np.nanargmin, min: np.nanargmin, np.min: np.nanargmin,
    'max': np.nanargmax, max: np.nanargmax, np.max: np.nanargmax
})

Mode = Literal['uniform', 'cluster', 'random']


def distribute_idx(core: ArrayLike, idx: Union[int, Iterable[int]], f: float,
                   mode: Mode = 'uniform', **kwargs: Any) -> np.ndarray:
    r"""Create a new distribution of atomic indices from **idx** of length :code:`f * len(idx)`.

    Parameters
    ----------
    core : array-like [:class:`float`], shape :math:`(m, 3)`
        A 2D array-like object (such as a :class:`Molecule` instance) consisting
        of Cartesian coordinates.

    idx : :class:`int` or :class:`Iterable<collections.abc.Iterable>` [:class:`int`], shape :math:`(i,)`
        An integer or iterable of unique integers representing the 0-based indices of
        all anchor atoms in **core**.

    f : :class:`float`
        A float obeying the following condition: :math:`0.0 < f \le 1.0`.
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
    :class:`numpy.ndarray` [:class:`int`], shape :math:`(f*i,)`
        A 1D array of atomic indices.
        If **idx** has :math:`i` elements,
        then the length of the returned list is equal to :math:`\max(1, f*i)`.

    See Also
    --------
    :func:`uniform_idx`
        Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    :func:`cluster_idx`
        Return the column-indices of **dist** which yield a clustered distribution.

    """  # noqa
    # Convert **idx** into an array
    idx_arr = as_nd_array(idx, dtype=int)

    # Validate the input
    if mode not in MODE_SET:
        raise ValueError(f"Invalid value for 'mode' ({reprlib.repr(mode)}); "
                         f"accepted values: {reprlib.repr(tuple(MODE_SET))}")
    elif not (0.0 < f <= 1.0):
        raise ValueError("'f' should be larger than 0.0 and smaller than or equal to 1.0; "
                         f"observed value: {reprlib.repr(f)}")

    # Create an array of indices
    stop = max(1, int(round(f * len(idx_arr))))
    if mode in ('uniform', 'cluster'):
        xyz = np.array(core, dtype=float, ndmin=2, copy=False)[idx_arr]
        dist = edge_dist(xyz) if kwargs.get('follow_edge', False) else cdist(xyz, xyz)
        operation = 'min' if mode == 'uniform' else 'max'
        generator1 = uniform_idx(dist, operation=operation,
                                 start=kwargs.get('start', None),
                                 cluster_size=kwargs.get('cluster_size', 1),
                                 randomness=kwargs.get('randomness', None),
                                 weight=kwargs.get('weight', lambda x: np.exp(-x)))
        generator2 = islice(generator1, stop)
        ret = idx_arr[np.fromiter(generator2, count=stop, dtype=int)]

    elif mode == 'random':
        ret = np.random.permutation(idx_arr)

    # Return a list of `p * len(idx_arr)` atomic indices
    return ret[:stop]


def uniform_idx(dist: np.ndarray, operation: str = 'min',
                cluster_size: Union[int, Iterable[int]] = 1,
                start: Optional[int] = None, randomness: Optional[float] = None,
                weight: Callable[[np.ndarray], np.ndarray] = lambda x: np.exp(-x)
                ) -> Generator[int, None, None]:
    r"""Yield the column-indices of **dist** which yield a uniform or clustered distribution.

    Given the (symmetric) distance matrix :math:`\boldsymbol{D} \in \mathbb{R}^{n,n}` and
    the vector :math:`\boldsymbol{a} \in \mathbb{N}^{m}`
    (representing a subset of indices in :math:`D`),
    then the :math:`i`'th element :math:`a_{i}` is
    defined below. All elements of :math:`\boldsymbol{a}` are
    furthermore constrained to be unique.
    :math:`f(x)` is herein a, as of yet unspecified,
    function for weighting each individual distance.

    Following the convention used in python, the :math:`\boldsymbol{X}[0:3, 1:5]` notation is
    herein used to denote the submatrix created
    by intersecting rows :math:`0` up to (but not including) :math:`3` and
    columns :math:`1` up to (but not including) :math:`5`.

    .. math::

        \DeclareMathOperator*{\argmin}{\arg\!\min}
        a_{i} = \begin{cases}
            \argmin\limits_{k \in \mathbb{N}} \sum f \bigl( \boldsymbol{D}_{k,:} \bigr) &
            \text{if} & i=0 \\
            \argmin\limits_{k \in \mathbb{N}}
                \sum f \bigl( \boldsymbol{D}[k, \boldsymbol{a}[0:i]] \bigr) &
            \text{if} & i > 0
        \end{cases}

    Default weighting function: :math:`f(x) = e^{-x}`.

    The row in :math:`D` corresponding to :math:`a_{0}`
    can alternatively be specified by **start**.

    The :math:`\text{argmin}` operation can be exchanged for :math:`\text{argmax}` by setting
    **operation** to ``"max"``, thus yielding a clustered- rather than uniform-distribution.

    The **cluster_size** parameter allows for the creation of uniformly
    distributed clusters of size :math:`r`.
    Herein the vector of indices, :math:`\boldsymbol{a} \in \mathbb{N}^{m}` is
    for the purpose of book keeping reshaped
    into the matrix :math:`\boldsymbol{A} \in \mathbb{N}^{q, r} \; \text{with} \; q*r = m`.
    All elements of :math:`\boldsymbol{A}` are, again, constrained to be unique.

    .. math::

        \DeclareMathOperator*{\argmin}{\arg\!\min}
        A_{i,j} = \begin{cases}
            \argmin\limits_{k \in \mathbb{N}} \sum f \bigl( \boldsymbol{D}_{k,:} \bigr) &
            \text{if} & i=0; \; j=0 \\
            \argmin\limits_{k \in \mathbb{N}}
                \sum f \bigl( \boldsymbol{D}[k; \boldsymbol{A}[0:i, 0:r] \bigl) &
            \text{if} & i > 0; \; j = 0 \\
            \argmin\limits_{k \in \mathbb{N}}
                \dfrac{\sum f \bigl( \boldsymbol{D}[k, \boldsymbol{A}[0:i, 0:r] \bigr)}
                {\sum f \bigl( \boldsymbol{D}[k, \boldsymbol{A}[i, 0:j] \bigr)} &
            \text{if} & j > 0
        \end{cases}

    Parameters
    ----------
    dist : :class:`numpy.ndarray` [:class:`float`], shape :math:`(n, n)`
        A symmetric 2D NumPy array (:math:`D_{i,j} = D_{j,i}`) representing the
        distance matrix :math:`D`.

    operation : :class:`str`
        Whether to minimize or maximize the distance between points.
        Accepted values are ``"min"`` and ``"max"``.

    start : :class:`int`, optional
        The index of the starting row in **dist**.
        If ``None``, start in whichever row contains the global minimum
        (:math:`\DeclareMathOperator*{\argmin}{\arg\!\min} \argmin\limits_{k \in \mathbb{N}} ||\boldsymbol{D}_{k, :}||_{p}`) or maximum
        (:math:`\DeclareMathOperator*{\argmax}{\arg\!\max} \argmax\limits_{k \in \mathbb{N}} ||\boldsymbol{D}_{k, :}||_{p}`).
        See **operation**.

    cluster_size : :class:`int` or :class:`Iterable<collections.abc.Iterable>` [:class:`int`]
        An integer or iterable of integers representing the size of clusters.
        Used in conjunction with :code:`operation = "max"` for creating a uniform distribution
        of clusters.
        :code:`cluster_size = 1` is equivalent to a normal uniform distribution.

        Providing **cluster_size** as an iterable of integers will create clusters
        of varying, user-specified, sizes.
        For example, :code:`cluster_size = range(1, 4)` will continuesly create clusters
        of sizes 1, 2 and 3.
        The iteration process is repeated until all atoms represented by **dist** are exhausted.

    randomness : :class:`float`, optional
        If not ``None``, represents the probability that a random index
        will be yielded rather than obeying **operation**.
        Should obey the following condition: :math:`0 \le randomness \le 1`.

    weight : :data:`Callable<typing.Callable>`
        A callable for applying weights to the distance; default: :math:`e^{-x}`.
        The callable should take an array as argument and return a new array,
        *e.g.* :func:`numpy.exp`.

    Yields
    ------
    :class:`int`
        Yield the column-indices specified in :math:`\boldsymbol{d}`.

    """  # noqa
    # Truncate and square the distance matrix
    dist_sqr = np.array(dist, dtype=float, copy=True)
    np.fill_diagonal(dist_sqr, np.nan)
    dist_sqr = weight(dist_sqr)

    # Use either argmin or argmax
    try:
        arg_func = OPERATION_MAPPING[operation]
    except KeyError as ex:
        raise ValueError(f"Invalid value for 'operation' ({reprlib.repr(operation)}); "
                         "accepted values: ('min', 'max')") from ex
    start = arg_func(np.nansum(dist_sqr, axis=1)) if start is None else start

    if randomness is not None:
        arg_func = _parse_randomness(randomness, arg_func, len(dist))

    # Yield the first index
    try:
        dist_1d_sqr = dist_sqr[start].copy()
    except IndexError as ex:
        if not hasattr(start, '__index__'):
            raise TypeError("'start' expected an integer or 'None'; observed type: "
                            f"'{start.__class__.__name__}'") from ex
        raise ValueError("index 'start={start}' is out of bounds: 'len(dist)={len(dist)}'") from ex

    dist_1d_sqr[start] = np.nan
    yield start

    # Return a generator for yielding the remaining indices
    if cluster_size == 1:
        generator = _min_or_max(dist_sqr, dist_1d_sqr, arg_func)
    else:
        generator = _min_and_max(dist_sqr, dist_1d_sqr, arg_func, cluster_size)

    for i in generator:
        yield i


def _min_or_max(dist_sqr: np.ndarray, dist_1d_sqr: np.ndarray,
                arg_func: Callable[[np.ndarray], int]) -> Generator[int, None, None]:
    """Helper function for :func:`uniform_idx` if :code:`cluster_size == 1`."""
    for _ in range(len(dist_1d_sqr)-1):
        i = arg_func(dist_1d_sqr)
        dist_1d_sqr[i] = np.nan
        dist_1d_sqr += dist_sqr[i]
        yield i


def _min_and_max(dist_sqr: np.ndarray, dist_1d_sqr: np.ndarray,
                 arg_func: Callable[[np.ndarray], int],
                 cluster_size: Union[int, Iterable[int]] = 1) -> Generator[int, None, None]:
    """Helper function for :func:`uniform_idx` if :code:`cluster_size != 1`."""
    # Construct a boolean array denoting the start of new clusters
    bool_ar = np.zeros(len(dist_1d_sqr), dtype=bool)
    if isinstance(cluster_size, abc.Iterable):
        bool_indices = _parse_cluster_size(len(bool_ar), cluster_size)
        bool_ar[bool_indices] = True
    else:
        try:
            bool_ar[::cluster_size] = True
        except ValueError as ex:
            raise ValueError("'cluster_size' cannot be zero; oberved value: "
                             f"{reprlib.repr(cluster_size)}") from ex
        except TypeError as ex:
            raise TypeError("'cluster_size' expected a non-zero integer or iterable of integers; "
                            f"observed type: '{cluster_size.__class__.__name__}'") from ex
    bool_ar = bool_ar[1:]  # Skip the first element as it was already yielded in uniform_idx()

    dist_cluster = dist_1d_sqr.copy()
    for bool_ in bool_ar:
        if bool_:  # The start of a new cluster
            dist_1d_sqr += dist_cluster
            dist_cluster[:] = 0
            dist_1d = dist_1d_sqr
        else:  # The growth of a cluster
            dist_1d = dist_1d_sqr / dist_cluster

        j = arg_func(dist_1d)
        dist_1d_sqr[j] = np.nan
        dist_cluster += dist_sqr[j]
        yield j


def _random_arg_func(dist_1d: np.ndarray, arg_func: Callable[[np.ndarray], int],
                     threshold: float, idx: np.ndarray,
                     rand_func: Callable[[], float] = np.random.sample) -> int:
    """Return a random element from **idx** if :code:`threshold > rand_func()`, otherwise call :code:`arg_func(dist_1d)`.

    Elements in **dist_1d** which are `nan` will be used for masking **idx**.

    """  # noqa
    if threshold > rand_func():
        return np.random.choice(idx[~np.isnan(dist_1d)])
    return arg_func(dist_1d)


def _parse_cluster_size(ar_size: int, clusters: Iterable[int]) -> np.ndarray:
    """Return indices for all ``True`` values in the boolean array of :func:`_min_and_max`."""
    generator = takewhile(lambda x: x < ar_size, cycle_accumulate(clusters))
    try:
        return np.fromiter(generator, dtype=int)
    except TypeError as ex:
        raise TypeError("'cluster_size' expected a non-zero integer or iterable of integers; "
                        f"{ex}") from ex


def _parse_randomness(randomness: float, arg_func: Callable[[np.ndarray], int],
                      n: int) -> Callable[[np.ndarray], int]:
    """Modifiy **arg_func** such that there is a **randomness** chance to return a random index from the range **n**."""  # noqa
    try:
        assert (0 <= randomness <= 1)
    except TypeError as ex:
        raise TypeError("'randomness' expected a float larger than 0.0 and smaller than 1.0; "
                        f"observed type: '{randomness.__class__.__name__}'") from ex
    except AssertionError as ex:
        raise ValueError("'randomness' expected a float larger than 0.0 and smaller than 1.0; "
                         f"observed value: {reprlib.repr(randomness)}") from ex

    return functools.partial(_random_arg_func,
                             arg_func=arg_func,
                             threshold=randomness,
                             idx=np.arange(n))
