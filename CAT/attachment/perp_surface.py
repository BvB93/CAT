"""A module for constructing vectors perpendicular to a molecules surface.

Index
-----
.. currentmodule:: CAT.attachment.perp_surface
.. autosummary::
    get_surface_vec
    plot_vectors

API
---
.. autofunction:: get_surface_vec
.. autofunction:: plot_vectors

"""

from typing import Union, Optional, Any

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

from scm.plams import Molecule

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    PLT: Optional[ImportError] = None
    Figure = plt.Figure
except ImportError as ex:
    PLT = ex
    Figure = 'matplotlib.pyplot.Figure'

__all__ = ['get_surface_vec', 'plot_vectors']


def get_surface_vec(mol: Union[Molecule, np.ndarray],
                    anchor: Union[None, Molecule, np.ndarray] = None) -> np.ndarray:
    """Construct a set of vectors perpendicular to the surface of **mol** and assign one to each atom in **anchor**.

    Utilizes a convex hull algorithm for identifying and partitioning the surface.

    Parameters
    ----------
    mol : array-like [:class:`float`], shape :math:`(n, 3)`
        A 2D array-like object with the Cartesian coordinates of a molecule.

    anchor : array-like [:class:`float`], shape :math:`(m, 3)`, optional
        A 2D array-like object with the Cartesian coordinates of a set of anchor atoms.
        If ``None``, default to **mol**.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`float`], shape :math:`(m, 3)`
        An array with vectors perpendicular to the surface of **mol**,
        one for each atom in **anchor**.

    See Also
    --------
    :class:`ConvexHull<scipy.spatial.ConvexHull>`
        Convex hulls in N dimensions.

    """  # noqa
    xyz = np.array(mol, dtype=float, ndmin=2, copy=False)
    if anchor is None:
        anchor = xyz
    else:
        anchor = np.array(anchor, dtype=float, ndmin=2, copy=False)

    # Construct the convex hull and extract the vertices
    hull = ConvexHull(xyz)
    simplice = np.swapaxes(xyz[hull.simplices], 0, 1)
    simplice_center = simplice.mean(axis=0)

    # Construct and return the surface vectors
    vec = _get_perp_vecs(*simplice)
    _flip_vec(simplice_center, vec)
    return vec[_find_nearest_center(anchor, simplice_center)]


def plot_vectors(vec: np.ndarray,
                 xyz: Optional[np.ndarray] = None,
                 show: bool = True, **kwargs: Any) -> Figure:
    r"""Create a 3D plot of all (3D) vectors in **vec**.

    Parameters
    ----------
    vec : array-like [:class:`float`], shape :math:`(n, 3)`
        A 2D array-like object representing :math:`n` vectors.

    xyz : array-like [:class:`float`], shape :math:`(n, 3)`, optional
        An array with the Cartesian coordinates defining the
        origin of each vector in **vec**.
        If ``None``, default to the origin (:code:`[0, 0, 0]`).

    show : :class:`bool`
        Show the created figure.

    \**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for
        :meth:`Axes.quiver()<matplotlib.pyplot.Axes.quiver>`
        such as the **length** keyword.

    Returns
    -------
    :class:`Figure<matplotlib.pyplot.Figure>`
        The resulting matplotlib Figure.

    """
    if PLT is not None:
        raise PLT

    # Parse arguments
    vec = np.array(vec, ndmin=2, dtype=float, copy=False)
    if xyz is not None:
        xyz = np.array(xyz, ndmin=2, dtype=float, copy=False)
    else:
        xyz = np.array(0.0, ndmin=2)

    # Extract the x, y and z coordinates
    if xyz.size == 1:
        x = y = z = xyz
    else:
        x, y, z = xyz.T

    # Extract the x, y and z components of the vector
    u, v, w = vec.T

    # Construct the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    _set_axis_limit(ax, x, y, z, u, v, w)
    ax.quiver(x, y, z, u, v, w, **kwargs)

    if show:
        plt.show(block=True)
    return fig


def _set_axis_limit(ax, x, y, z, u, v, w) -> None:
    """Set the axis limits for :func:`plot_vectors`."""
    x_mean = x.mean()
    y_mean = y.mean()
    z_mean = z.mean()
    dx = max(1, max(max(abs(x)), max(abs(u))) - x_mean)
    dy = max(1, max(max(abs(y)), max(abs(v))) - y_mean)
    dz = max(1, max(max(abs(z)), max(abs(w))) - z_mean)
    ax.set_xlim([x_mean-dx, x_mean+dx])
    ax.set_ylim([y_mean-dy, y_mean+dy])
    ax.set_zlim([z_mean-dz, z_mean+dz])


def _find_nearest_center(anchor: np.ndarray, center: np.ndarray) -> np.ndarray:
    """Find the points in **center** closest to those in **anchor**."""
    tree = cKDTree(center)
    _, idx = tree.query(anchor, k=1)
    return idx


def _get_perp_vecs(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Construct a unit vector perpendicular to a set of triangular polygons."""
    v1 = p2 - p1
    v2 = p3 - p1
    vec = np.cross(v1, v2)
    vec /= np.linalg.norm(vec, axis=-1)[..., None]
    return vec


def _flip_vec(simplice_centra: np.ndarray, perp_vec: np.ndarray) -> None:
    """Ensure that **perp_vec** is pointing away from the surface defined by **simplice_centra**."""
    center = simplice_centra.mean(axis=0)
    dot = np.einsum('ij,ij->i', simplice_centra - center, perp_vec)
    perp_vec[dot < 0] *= -1
