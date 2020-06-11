"""Functions for creating distance matrices.

Index
-----
.. currentmodule:: CAT.attachment.edge_distance
.. autosummary::
    edge_dist
    array_combinations
    to_convex

API
---
.. autofunction:: edge_dist
.. autofunction:: array_combinations
.. autofunction:: to_convex

"""

import reprlib
from typing import Optional, Any

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import ConvexHull

from nanoutils import array_combinations

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    PLT: Optional[ImportError] = None
    Figure = plt.Figure
except ImportError as ex:
    PLT = ex
    Figure = 'matplotlib.pyplot.Figure'

__all__ = ['edge_dist', 'plot_polyhedron']


def to_convex(xyz: np.ndarray, n: float = 1.0) -> np.ndarray:
    r"""Round all edges in **xyz** by a factor **n**: (:math:`0 < n \le 1`)."""
    if not (0.0 < n <= 1.0):
        raise ValueError(f"Expected '0 < n <= 1'; observed: {reprlib.repr(n)}")
    xyz = np.asarray(xyz, dtype=float)
    xyz = xyz - xyz.mean(axis=0)

    r = np.linalg.norm(xyz, axis=-1)
    dr = n * (r.mean() - r)

    elongate = (dr + r) / r
    xyz *= elongate[..., None]
    return xyz


def edge_dist(xyz: np.ndarray, n: float = 1.0,
              edges: Optional[np.ndarray] = None) -> np.ndarray:
    r"""Calculate all shortest paths between all points in the polyhedron **xyz** by traversing its edges.

    After converting **xyz** into a polyhedron with triangular faces,
    using a convex hull algorithm,
    the shortest paths between all possible point pairs is calculated using
    the Dijkstra algorithm.
    Paths are constrained to the edges of the polyhedron,
    as opposed to a "normal" distance matrix each element thus represents
    the shortest path along a set of 1D lines rather than a 3D volume.

    Given the matrix of Cartesian coordinates :math:`\boldsymbol{X} \in \mathbb{R}^{n, 3}`,
    the matching distance matrix :math:`\boldsymbol{D} \in \mathbb{R}^{n, n}` and
    the matching edge-distance matrix :math:`\boldsymbol{D}^{\text{edge}} \in \mathbb{R}^{n, n}`,
    then element :math:`D_{i,j}` is defined as following:

    .. math::

        D_{i, j} = ||\boldsymbol{X}_{i,:} - \boldsymbol{X}_{j,:}||

    while :math:`D_{i, j}^{\text{edge}}` is defined as:

    .. math::

        D_{i, j}^{\text{edge}} = \min_{\boldsymbol{p} \in \mathbb{N}^{m}; m \in \mathbb{N}}
        \sum_{k=0}^{m-1} || X_{p_{k},:} - X_{p_{k+1},:} ||
        \quad \text{with} \quad p_{0} = i \quad \text{and} \quad p_{m} = j

    The vector :math:`\boldsymbol{p}` is a path,
    represented by the indices of neighbouring vertices in :math:`\boldsymbol{X}`,
    whose length (*i.e.* sum of Euclidean distances) is to-be minimized.

    Notes
    -----
    All points in **xyz** are projected on the surface of a sphere during the construction
    of the convex hull (if :code:`n != 0.0`);
    the quality of the constructed polyhedron will thus depend on the "convexness" of **xyz**.
    For highly concave structures (*e.g.* a torus) it is strongly recommended to
    manually pass all edge indices using the **edges** parameter.

    Parameters
    ----------
    xyz : array-like [:class:`float`], shape :math:`(m, 3)`
        A 2D array-like object of Cartesian coordinates representing a polyhedron.
        The supplied polyhedron should be convex in shape.

    n : :class:`float`
        Smoothing factor for constructing a convex hull.
        Should obey :math:`0 <= n <= 1`.

    edges : array-like [:class:`int`], shape :math:`(n, 2)`, optional
        A 2D array-like object with all indice-pairs in **xyz** representing polyhedron edges.

    Returns
    -------
    :math:`(m, m)` :class:`numpy.ndarray`
        A 2D array containing all possible (Euclidean) distance-pairs in **xyz**.
        Distances are calculated by traversing the shortest path along the edges of **xyz**,
        rather than moving directly through space.

    See Also
    --------
    :class:`ConvexHull<scipy.spatial.ConvexHull>`
        Convex hulls in N dimensions.

    :func:`dijkstra<scipy.sparse.csgraph.dijkstra>`
        Dijkstra algorithm using Fibonacci Heaps.

    :func:`cdist<scipy.spatial.distance.cdist>`
        Compute distance between each pair of the two collections of inputs.

    """  # noqa
    xyz = np.array(xyz, dtype=float, ndmin=2, copy=False)
    xyz_len = len(xyz)

    # Create an array with all index-pairs forming the polyhedron edges
    if edges is None:
        xyz_convex = to_convex(xyz, n=n) if not np.allclose(n, 0.0) else xyz
        hull = ConvexHull(xyz_convex)
        idx = array_combinations(hull.simplices, r=2)
        idx.shape = -1, 2
    else:
        idx = np.asarray(edges, dtype=int)
    i, j = idx.T

    # Create a distancea matrix containing all edge distances
    dist = np.zeros((xyz_len, xyz_len), dtype=float)
    dist[i, j] = dist[j, i] = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    dist_sparse = csr_matrix(dist)

    # Traverse the edges
    return dijkstra(dist_sparse,
                    directed=False,
                    indices=np.arange(xyz_len),
                    return_predecessors=False)


def plot_polyhedron(xyz: np.ndarray, triangles: Optional[np.ndarray] = None,
                    show: bool = True, **kwargs: Any) -> Figure:
    r"""Plot a polyhedron, represented by an array of Cartesian coordinates, with matplotlib.

    Parameters
    ----------
    xyz : array-like [:class:`float`], shape :math:`(m, 3)`
        A 2D array-like object representing the Cartesian coordinates of a polyhedron.

    triangles : array-like [:class:`int`], shape :math:`(n, 3)`, optional
        A 2D array-like object with all indice-pairs in **xyz** representing the triangular
        faces of the **xyz** polygon.

    show : :class:`bool`
        Show the created figure.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for
        :meth:`Axes.plot_trisurf()<matplotlib.pyplot.Axes.plot_trisurf>`.

    Returns
    -------
    :class:`Figure<matplotlib.pyplot.Figure>`
        The resulting matplotlib Figure.

    """
    if PLT is not None:
        raise PLT

    if 'cmap' not in kwargs:
        kwargs['cmap'] = plt.cm.Spectral

    xyz = np.asarray(xyz, dtype=float)
    triangles = ConvexHull(xyz).simplices if triangles is None else np.asarray(triangles, dtype=int)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf(*xyz.T, triangles=triangles, **kwargs)

    if show:
        plt.show(block=True)
    return fig
