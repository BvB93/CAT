"""A module designed for attaching ligands to cores.

Index
-----
.. currentmodule:: CAT.attachment.ligand_attach
.. autosummary::
    init_qd_construction
    construct_mol_series
    _read_database
    _get_indices
    _get_df
    ligand_to_qd
    _get_rotmat1
    _get_rotmat2
    rot_mol
    rot_mol_angle
    array_to_qd
    sanitize_dim_2
    sanitize_dim_3

API
---
.. autofunction:: init_qd_construction
.. autofunction:: construct_mol_series
.. autofunction:: _read_database
.. autofunction:: _get_indices
.. autofunction:: _get_df
.. autofunction:: ligand_to_qd
.. autofunction:: _get_rotmat1
.. autofunction:: _get_rotmat2
.. autofunction:: rot_mol
.. autofunction:: rot_mol_angle
.. autofunction:: array_to_qd
.. autofunction:: sanitize_dim_2
.. autofunction:: sanitize_dim_3

"""

from typing import List, Tuple, Any, Optional, NoReturn, Union, Iterable
from collections import abc

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from scm.plams import Molecule, Atom, Settings, MoleculeError
from assertionlib.ndrepr import aNDRepr

from .perp_surface import get_surface_vec
from ..mol_utils import get_index, round_coords  # noqa: F401
from ..workflows import WorkFlow, HDF5_INDEX, MOL, OPT, LIGAND_COUNT
from ..settings_dataframe import SettingsDataFrame
from ..data_handling import mol_to_file, WARN_MAP

__all__ = ['init_qd_construction']


def init_qd_construction(ligand_df: SettingsDataFrame, core_df: SettingsDataFrame,
                         construct_qd: bool = True) -> SettingsDataFrame:
    """Initialize the quantum dot construction.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of ligands.

    core_df : |CAT.SettingsDataFrame|_
        A dataframe of cores.

    construct_qd : :class:`bool`
        If ``False``, only construct and return the dataframe without filling it with molecules.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of quantum dots.

    """
    qd_df = _get_df(core_df.index, ligand_df.index, ligand_df.settings)
    qd_df[MOL] = None
    qd_df[LIGAND_COUNT] = 0
    qd_df.sort_index(inplace=True)
    if not construct_qd:
        return qd_df

    workflow = WorkFlow.from_template(qd_df, name='qd_attach')
    workflow.keep_files = False

    # Pull from the database; push unoptimized structures
    df_bool = workflow.from_db(qd_df, LIGAND_COUNT[0], read_mol=True)

    # Start the ligand optimization
    idx = df_bool[MOL]
    workflow(construct_mol_series, qd_df, columns=MOL, index=idx,
             core_df=core_df, ligand_df=ligand_df)

    for k, mol in core_df[MOL].items():
        qd_df.loc[k, LIGAND_COUNT] = len(mol.properties.dummies)
    workflow.to_db(qd_df, df_bool, columns=[MOL, LIGAND_COUNT], status='no_opt')

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    mol_format = qd_df.settings.optional.database.mol_format
    if mol_format and not qd_df.settings.optional.qd.optimize:
        path = workflow.path
        mol_array = qd_df.loc[idx, MOL].values
        mol_to_file(mol_array, path, mol_format=mol_format)

    return qd_df


def construct_mol_series(qd_df: SettingsDataFrame, core_df: pd.DataFrame,
                         ligand_df: pd.DataFrame, path: str,
                         allignment: str = 'sphere', **kwargs: Any) -> pd.Series:
    """Construct a Series of new quantum dots."""
    def _get_mol(i, j, k, m):
        ij = i, j
        km = k, m
        return ligand_to_qd(core_df.at[ij, MOL], ligand_df.at[km, MOL], path, allignment=allignment)

    mol_list = [_get_mol(i, j, k, m) for i, j, k, m in qd_df.index]
    return pd.Series(mol_list, index=qd_df.index, name=MOL, dtype=object)


def _get_indices(mol: Molecule,
                 index: Tuple[str, str, str, str]) -> List[int]:
    """Return a list with the indices of all atoms in the core plus ligand anchor atoms.

    Ligand anchor atoms are furthermore marked with the properties.anchor attribute.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    index : |tuple|_ [|str|_]
        A tuple of 4 strings.

    Returns
    -------
    |list|_ [|int|_]
        A list of atomic indices.

    """
    # Collect the indices of the atoms in the core
    ret = []
    for i, at in enumerate(mol, 1):
        if at.properties.pdb_info.ResidueName == 'COR':
            ret.append(i)
        else:
            break

    # Extract the index (within the ligand) of the ligand anchor atom
    index = index[3]
    for j, _ in enumerate(index):
        try:
            k = index[j:] - 1
            break
        except ValueError:
            pass
    k += i - 1

    # Append and return
    ref_name = mol[k+1].properties.pdb_info.Name
    for i, at in enumerate(mol.atoms[k:], k+1):
        if at.properties.pdb_info.Name == ref_name:
            at.properties.anchor = True
            ret.append(i)
    return ret


def _get_df(core_index: pd.MultiIndex,
            ligand_index: pd.MultiIndex,
            settings: Settings) -> SettingsDataFrame:
    """Create and return a new quantum dot dataframe.

    Parameters
    ----------
    core_index : |pd.MultiIndex|_
        A multiindex of the cores.

    ligand_index : |pd.MultiIndex|_
        A multiindex of the ligands.

    settings : |plams.Settings|_
        A Settings intance extracted from the ligand or core dataframe.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        An empty dataframe (*i.e.* filled with ``-1`` and ``False``) of quantum dots.

    """
    # Create the index
    idx_tups = [(i, j, k, l) for i, j in core_index for k, l in ligand_index]
    index = pd.MultiIndex.from_tuples(
        idx_tups, names=['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    )

    # Create the collumns
    column_tups = [HDF5_INDEX, OPT]
    columns = pd.MultiIndex.from_tuples(column_tups, names=['index', 'sub index'])

    # Create and return the quantum dot dataframe
    data = {HDF5_INDEX: 0, OPT: False}
    return SettingsDataFrame(data, index=index, columns=columns, settings=settings)


def ligand_to_qd(core: Molecule, ligand: Molecule, path: str,
                 allignment: str = 'sphere',
                 idx_subset: Optional[Iterable[int]] = None) -> Molecule:
    """Function that handles quantum dot (qd, *i.e.* core + all ligands) operations.

    Combine the core and ligands and assign properties to the quantom dot.

    Parameters
    ----------
    core : |plams.Molecule|_
        A core molecule.

    ligand : |plams.Molecule|_
        A ligand molecule.

    allignment : :class:`str`
        How the core vector(s) should be defined.
        Accepted values are ``"sphere"`` and ``"surface"``:

        * ``"sphere"``: Vectors from the core anchor atoms to the center of the core.
        * ``"surface"``: Vectors perpendicular to the surface of the core.

        Note that for a perfect sphere both approaches are equivalent.

    idx_subset : :class:`Iterable<collections.anc.Iterable>` [:class:`int`], optional
        An iterable with the (0-based) indices defining a subset of atoms in **core**.
        Only relevant in the construction of the convex hull when ``allignment=surface``.

    Returns
    -------
    |plams.Molecule|_
        A quantum dot consisting of a core molecule and *n* ligands

    """
    def get_name() -> str:
        core_name = core.properties.name
        anchor = str(qd[-1].properties.pdb_info.ResidueNumber - 1)
        lig_name = ligand.properties.name
        return f'{core_name}__{anchor}_{lig_name}'

    idx_subset = idx_subset if idx_subset is not None else ...

    # Define vectors and indices used for rotation and translation the ligands
    vec1 = np.array([-1, 0, 0], dtype=float)  # All ligands are already alligned along the X-axis
    idx = ligand.get_index(ligand.properties.dummies) - 1
    ligand.properties.dummies.properties.anchor = True

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    if allignment == 'sphere':
        vec2 = np.array(core.get_center_of_mass()) - sanitize_dim_2(core.properties.dummies)
        vec2 /= np.linalg.norm(vec2, axis=1)[..., None]
    elif allignment == 'surface':
        vec2 = -get_surface_vec(np.array(core)[idx_subset],
                                core.as_array(core.properties.dummies))
    else:
        raise ValueError(repr(allignment))

    lig_array = rot_mol(ligand, vec1, vec2, atoms_other=core.properties.dummies, core=core, idx=idx)
    qd = core.copy()
    array_to_qd(ligand, lig_array, mol_out=qd)
    qd.round_coords()

    # Set properties
    qd.properties = Settings({
        'indices': [i for i, at in enumerate(qd, 1) if
                    at.properties.pdb_info.ResidueName == 'COR' or at.properties.anchor],
        'path': path,
        'name': get_name(),
        'job_path': [],
        'prm': ligand.properties.get('prm')
    })

    # Print and return
    _evaluate_distance(qd, qd.properties.name)
    return qd


def _get_rotmat1(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """Calculate the rotation matrix for rotating **vec1** to **vec2**.

    Returns a unit matrix if the length of one of the vectors is 0.

    Paramaters
    ----------
    vec1 : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.
        Represents the initial orientation.

    vec2 : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.
        Represents the final orientation.

    Returns
    -------
    :math:`n*3*3` |np.ndarray|_
        An array consisting of :math:`n` rotation matrices.

    """
    # Increase the vector dimensionality if required
    u_vec1 = np.array(vec1, dtype=float, ndmin=2, copy=False)
    u_vec2 = np.array(vec2, dtype=float, ndmin=2, copy=False)

    with np.errstate(divide='raise', invalid='raise'):
        # Return a unit matrix if either vec1 or vec2 is zero
        try:
            u_vec1 /= np.linalg.norm(u_vec1, axis=1)[:, None]
            u_vec2 /= np.linalg.norm(u_vec2, axis=1)[:, None]
        except FloatingPointError:  # i.e. division by zero
            shape = max(len(u_vec1), len(u_vec2)), 3, 3
            ret = np.zeros(shape, dtype=float)
            ret[:] = np.identity(3)
            return ret

    v1, v2, v3 = np.cross(u_vec1, u_vec2).T
    v0 = np.zeros(max(len(u_vec1), len(u_vec2)))
    m = np.array([[v0, -v3, v2],
                  [v3, v0, -v1],
                  [-v2, v1, v0]]).T

    with np.errstate(invalid='ignore'):
        if len(u_vec1) > 1 and len(u_vec2) > 1:
            ret = np.identity(3) + m + ((m@m).T / (1 + u_vec2.T@u_vec1).T[..., None]).T
        elif len(u_vec1) == 1:
            ret = np.identity(3) + m + ((m@m).T / (1 + u_vec2@u_vec1[0].T).T).T
        elif len(u_vec2) == 1:
            ret = np.identity(3) + m + ((m@m).T / (1 + u_vec2[0]@u_vec1.T).T).T
        else:
            raise ValueError('vec1 and vec2 expect 1- or 2-dimensional array-like objects; '
                             f'observed shapes: vec1: {np.asarray(vec1).shape} and vec2: '
                             f'{np.asarray(vec2).shape}')

    # Ensure that there are no nan values in the rotation matrix
    isnan = np.isnan(ret)
    if not isnan.any():
        return ret
    else:
        i = np.unique(np.nonzero(isnan)[0])
        ret[i] = -np.eye(3)
        return ret


def _get_rotmat2(vec: np.ndarray, step: float = (1/16)) -> np.ndarray:
    r"""Calculate the rotation matrix for rotating m vectors along their axes, each vector yielding :math:`k = (2 / step)` possible rotations.

    Paramaters
    ----------
    vec : :math:`3` or :math:`n*3` |np.ndarray|_
        An array-like object representing a single or :math:`n` vectors of length 3.

    step : float
        The rotation stepsize as fraction of :math:`2*/pi`.

    """  # noqa
    # Increase the vector dimensionality if required and create unit vectors
    v = np.array(vec, dtype=float, ndmin=2, copy=False)
    v /= np.linalg.norm(v, axis=1)[:, None]

    v1, v2, v3 = v.T
    zero = np.zeros(len(v))
    w = np.array([[zero, -v3, v2],
                  [v3, zero, -v1],
                  [-v2, v1, zero]]).T

    step_range = np.pi * np.arange(0.0, 2.0, step)
    a1 = np.sin(step_range)[:, None, None, None]
    a2 = (1 - np.cos(step_range))[:, None, None, None]

    return np.identity(3) + a1 * w + a2 * w@w


def rot_mol(xyz_array: np.ndarray,
            vec1: np.ndarray,
            vec2: np.ndarray,
            idx: int = 0,
            atoms_other: Optional[np.ndarray] = None,
            core: Optional[np.ndarray] = None,
            bond_length: Optional[int] = None,
            step: float = 1/16,
            dist_to_self: bool = True,
            ret_min_dist: bool = False) -> np.ndarray:
    r"""Rotate **xyz_array**.

    Paramaters
    ----------
    xyz_array : :math:`m*n*3` |np.ndarray|_ or :math:`m` |list|_ [|plams.Molecule|_]:
        A 3D array-like object or list of :math:`m` PLAMS molecules consisting of :math:`n` atoms.
        This array represents the to-be rotated molecules.

    vec1 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the initial orientation of **xyz_array**.

    vec2 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the final orientation of **xyz_array**.
        Also used as a rotation axis.

    idx : int:
        An atomic index or sequence of atomic indices.
        Used for translational purposes (see **atoms_other** and **bond_length**).
        Expects 0-based indices.

    atoms_other : :math:`n*3` |np.ndarray|_:
        Optional: A 2D array of Cartesian coordinates.
        If not ``None`` then **xyz_array** will be translated to ensure
        :code:`xyz_array[:, idx] == atoms_other`.

    bond_length : float
        Optional: One or multiple floats representing bond lengths.
        Performs a translation of :code:`xyz_array[:, idx]` along **vec2** by **bond_length**.

    step : float
        The rotation stepsize as fraction of :math:`2*\pi`.

    dist_to_self : bool
        pass

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_:
        The rotated **xyz_array**.

    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    xyz = sanitize_dim_3(xyz_array)
    vec1 = sanitize_dim_2(vec1)
    vec2 = sanitize_dim_2(vec2)

    # Define slices
    idx1 = np.arange(len(xyz)), idx
    idx2 = np.arange(len(xyz)), slice(None), idx

    # Translate xyz[idx] to the origin and rotate
    xyz -= xyz[idx1][..., None, :]
    rotmat1 = _get_rotmat1(vec1, vec2)
    xyz = xyz@rotmat1

    # Create all k possible rotations of all m ligands
    rotmat2 = _get_rotmat2(vec2, step=step)
    xyz = np.swapaxes(xyz@rotmat2, 0, 1)
    if atoms_other is None:
        return xyz

    # Translate the the molecules in xyz_array
    at_other = sanitize_dim_2(atoms_other)
    xyz += (at_other[..., None, :] - xyz[idx2])[..., None, :]
    if bond_length is not None:
        bond_length = np.asarray(bond_length)
        mult = (bond_length / np.linalg.norm(vec2, axis=1))[:, None]
        xyz -= (vec2 * mult)[:, None, None, :]

    # Returns the conformation of each molecule that maximizes the inter-moleculair distance
    # Or return all conformations if dist_to_self = False and atoms_other = None
    return rotation_check_kdtree(xyz, at_other, core, ret_min_dist=ret_min_dist)


def rotation_check_kdtree(xyz: np.ndarray, core_anchor: np.ndarray, core: np.ndarray,
                          k: int = 10,
                          ret_min_dist: bool = False):
    """Perform the rotation check using SciPy's :class:`cKDTree<scipy.spatial.cKDTree.

    Parameters
    ----------
    xyz : :math:`(m, n, l, 3)` :class:`numpy.ndarray`
        A 4D array of Cartesian coordinates representing the :math:`n` rotameters of :math:m`
        ligands with :math:`l` atoms each.

    at_other : :math:`(m, 3)` :class:`numpy.ndarray`
        A 2D array with the Cartesian of neighbouring atoms.

    k : :class:`int`
        The number of to-be considered neighbours when performing the rotation check.

    See Also
    --------
    :meth:`cKDTree.query<scipy.spatial.cKDTree.query>`
        Query the kd-tree for nearest neighbors.

    """
    a, b, c, d = xyz.shape
    ret = np.empty((a, c, d), order='F')
    min_dist = np.empty(len(ret))
    distance_upper_bound = _get_distance_upper_bound(core_anchor)

    # Shrink down the core, keep the 6 atoms closest to the core-anchor
    core = np.asarray(core)
    tree = cKDTree(core)
    _, idx = tree.query(core_anchor, k=range(2, 7), distance_upper_bound=10)
    core_subset = core[np.unique(idx.ravel())]

    at_other = core_subset
    for i, ar in enumerate(xyz):
        tree = cKDTree(at_other)
        dist, _ = tree.query(ar.reshape(b*c, d), k=k, distance_upper_bound=distance_upper_bound)
        dist.shape = b, c, k

        weighted_dist = np.exp(-dist).sum(axis=(1, 2))
        idx_min = weighted_dist.argmin()
        at_other = np.concatenate((at_other, ar[idx_min]))

        ret[i] = ar[idx_min]
        min_dist[i] = weighted_dist[idx_min]

    if ret_min_dist:
        return ret, min_dist
    return ret


def _evaluate_distance(mol: np.ndarray, name: str,
                       threshold: float = 1.0,
                       action: str = 'warn') -> Union[None, NoReturn]:
    """Eavluate all the distance matrix of **xyz3D** and perform **action** when distances are below **threshold**."""  # noqa
    try:
        action_func = WARN_MAP[action]
    except KeyError as ex:
        raise ValueError("'action' expected either 'warn', 'raise' or 'ignore'; "
                         f"observed value: {action!r}") from ex
    else:
        if action == 'ignore':
            return None

    xyz = np.asarray(mol)

    tree = cKDTree(xyz)
    dist, idx = tree.query(xyz, k=[2])
    dist = dist.T[0]
    idx = idx.T[0]

    bool_ar = dist < threshold
    if bool_ar.any():
        _idx2 = np.stack([np.arange(len(idx)), idx]).T
        _idx2 += 1
        _idx2.sort(axis=1)

        idx2 = np.unique(_idx2[bool_ar], axis=0)
        n = len(idx2)

        exc = MoleculeError(
            f"\nEncountered >= {n} unique atom-pairs at a distance shorter than "
            f"{threshold} Angstrom in {name!r}:\n{aNDRepr.repr(idx2.T)}"
        )
        action_func(exc)


def _get_distance_upper_bound(at_other: np.ndarray, r_min: float = 5.0,
                              r_max: float = 10.0) -> float:
    """Construct an estimate for **distance_upper_bound** based on the mean nearest-neighbour distance in **at_other**.

    The to-be returned value is clipped, if necessary, by **r_min** and **r_max**.

    """  # noqa
    tree = cKDTree(at_other)
    dist, _ = tree.query(at_other, k=2, distance_upper_bound=10.0)
    dist = dist.T[1]
    dist[dist == np.inf] = np.nan

    # Check if distance_upper_bound=10.0 is too small; try again with the full dist matrix if so
    if np.isnan(dist).all():
        dist = cdist(at_other, at_other)
        np.fill_diagonal(dist, np.inf)
        r = dist.min(axis=0).mean()
    else:
        r = np.nanmean(dist)

    # Clip and return
    return np.clip(r, r_min, r_max)


def rot_mol_angle(xyz_array: np.ndarray,
                  vec1: np.ndarray,
                  vec2: np.ndarray,
                  idx: int = 0,
                  atoms_other: Optional[np.ndarray] = None,
                  bond_length: Optional[float] = None) -> np.ndarray:
    """Rotate one or more molecules (**xyz_array**) from **vec1** to **vec2**.

    Paramaters
    ----------
    xyz_array : :math:`m*n*3` |np.ndarray|_ or :math:`m` |list|_ [|plams.Molecule|_]:
        A 3D array-like object or list of :math:`m` PLAMS molecules consisting of :math:`n` atoms.
        This array represents the to-be rotated molecules.

    vec1 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the initial orientation of **xyz_array**.

    vec2 : :math:`n*3` |np.ndarray|_:
        One or multiple vectors of length 3.
        Represents the final orientation of **xyz_array**.

    idx : int:
        An atomic index or sequence of atomic indices.
        Used for translational purposes (see **atoms_other** and **bond_length**).
        Expects 0-based indices.

    atoms_other : :math:`n*3` |np.ndarray|_:
        Optional: A 2D array of Cartesian coordinates.
        If not ``None`` then **xyz_array** will be translated to ensure
        :code:`xyz_array[:, idx] == atoms_other`.

    bond_length : float
        Optional: One or multiple floats representing bond lengths.
        Performs a translation of :code:`xyz_array[:, idx]` along **vec2** by **bond_length**.

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_:
        The rotated **xyz_array**.

    """
    # Turn all arguments into numpy arrays with appropiate dimensions
    _xyz_array = sanitize_dim_3(xyz_array)
    _vec1 = sanitize_dim_2(vec1)
    _vec2 = sanitize_dim_2(vec2)

    # Rotate and translate all n ligands; readjust bond lengths if bond_length is set
    rotmat = _get_rotmat1(_vec1, _vec2)
    ret = _xyz_array@rotmat
    _idx = slice(0, None), idx

    # Translate the the molecules in ret
    if atoms_other is not None:
        atoms_other = sanitize_dim_2(atoms_other)
        ret += (atoms_other - ret[_idx])[:, None, :]
        if bond_length is not None:
            mult = (np.asarray(bond_length) / np.linalg.norm(vec2, axis=1))[:, None]
            ret -= (vec2 * mult)[:, None, :]

    # Return a m*n*3 array
    return ret


def array_to_qd(mol: Molecule, xyz_array: np.ndarray,
                mol_out: Optional[Molecule] = None) -> Optional[List[Molecule]]:
    """Create :math:`n` copies of **mol** and update their Cartesian coordinates with **xyz_array**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A template PLAMS molecule consisting of :math:`n` atoms.

    xyz_array : :math:`m*n*3` |np.ndarray|_:
        A 3D array-like object representing the cartesian coordinates of
        :math:`m` molecules each with :math:`n` atoms.

    mol_out : |plams.Molecule|_
        Optional: If not ``None`` merge the to-be returned molecules with **mol_out**.

    Returns
    -------
    |list|_ [|Molecule|_]
        Optional: if **mol_out** = ``None`` return a list **mol** copies, each with its
        Cartesian coordinates replaced with a set of coordinates from **xyz_array**.

    """
    mol_list = []
    xyz_array = sanitize_dim_3(xyz_array)
    if mol_out is None:
        start = 2
    else:
        start = 1 + mol_out[-1].properties.get('pdb_info', {}).get('ResidueNumber', 1)

    for i, xyz in enumerate(xyz_array, start=start):
        mol_cp = mol.copy()
        mol_cp.from_array(xyz)
        for at in mol_cp:
            at.properties.pdb_info.ResidueNumber = i
        mol_list.append(mol_cp)

    # return a list of molecules or merge the list with an existing molecule
    if mol_out is None:
        return mol_list
    else:
        for m in mol_list:
            mol_out.add_molecule(m)
        return None


def _is_sequence(item) -> bool:
    return isinstance(item, abc.Sequence)


def _is_atom(item) -> bool:
    return isinstance(item, Atom)


def _is_mol(item) -> bool:
    return isinstance(item, Molecule)


def _is_atom_sequence(item) -> bool:
    return _is_mol(item) or (_is_sequence(item) and _is_atom(item[-1]))


def _is_mol_sequence(item) -> bool:
    return _is_sequence(item) and _is_atom_sequence(item[-1])


def sanitize_dim_2(value: Any) -> np.ndarray:
    """Convert a PLAMS atom or sequence of :math:`n` PLAMS atoms into a :math:`n*3` array.

    Parameters
    ----------
    value : object
        The to be parsed object.
        Acceptable types are:

        * A PLAMS :class:`Atom`.
        * A PLAMS :class:`Molecule`.
        * A Sequence consisting of PLAMS :class:`Atom`.
        * An array-like object with a dimensionality smaller than 2 and float-compatible elements.

    Returns
    -------
    :math:`n*3` |np.ndarray|_ [|np.float64|_]
        A 2D array consisting of floats.

    Raises
    ------
    ValueError
        Raised if dimensionality of the to-be returned array is higher than 2
        or the content of **arg** cannot be converted into an array of floats.

    """
    value = value.coords if isinstance(value, Atom) else value
    try:
        ret = np.array(value, dtype=float, copy=False)
    except TypeError:
        ret = Molecule.as_array(None, atom_subset=value)

    if ret.ndim < 2:
        ret.shape = (2 - ret.ndim) * (1,) + ret.shape
    elif ret.ndim > 2:
        raise ValueError(f"Failed to create a 2D array; observed dimensionality: {ret.ndim}")
    return ret


def sanitize_dim_3(value: Any, padding: float = np.nan) -> np.ndarray:
    """Convert a Molecule or sequence of :math:`m` molecules into an :math:`m*n*3` array.

    If necessary, the to-be returned array is padded with **padding** .

    Parameters
    ----------
    arg : object
        The to be parsed object.
        Acceptable types are:

        * A PLAMS :class:`Atom`
        * A PLAMS :class:`Molecule`
        * A (nested) Sequences consisting of PLAMS :class:`Atom`
        * A (nested) Sequences consisting of PLAMS :class:`Molecule`
        * An array-like object with a dimensionality smaller than 3 and float-compatible elements.

    padding : float
        A value used for padding the to-be returned array.
        Only relevant if **arg** consists of multiple molecule with different numbers of atoms.

    Returns
    -------
    :math:`m*n*3` |np.ndarray|_ [|np.float64|_]
        A 3D array consisting of floats.

    Raises
    ------
    ValueError
        Raised if dimensionality of the to-be returned array is higher than 3
        or the content of **arg** cannot be converted into an array of floats.

    """
    if _is_atom(value):
        return np.array(value.coords)[None, None, :]

    elif _is_atom_sequence(value):
        return Molecule.as_array(None, atom_subset=value)[None, :]

    elif _is_mol_sequence(value):
        max_at = max(len(mol) for mol in value)
        ret = np.full((len(value), max_at, 3), padding, order='F')
        for i, mol in enumerate(value):
            j = len(mol)
            ret[i, :j] = Molecule.as_array(None, atom_subset=mol)
        return ret

    else:
        ret = np.array(value, ndmin=3, dtype=float, copy=False)
        if ret.ndim > 3:
            raise ValueError(f"Failed to create a 3D array; observed dimensionality: {ret.ndim}")
        return ret
