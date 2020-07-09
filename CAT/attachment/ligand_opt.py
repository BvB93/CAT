"""A module designed for optimizing the geometry of ligands.

Index
-----
.. currentmodule:: CAT.attachment.ligand_opt
.. autosummary::
    init_ligand_opt
    _parse_overwrite
    read_data
    start_ligand_jobs
    optimize_ligand
    _ligand_to_db
    split_mol
    get_frag_size
    get_dihed
    set_dihed

API
---
.. autofunction:: init_ligand_opt
.. autofunction:: _parse_overwrite
.. autofunction:: read_data
.. autofunction:: start_ligand_jobs
.. autofunction:: optimize_ligand
.. autofunction:: _ligand_to_db
.. autofunction:: get_frag_size
.. autofunction:: split_bond
.. autofunction:: get_dihed
.. autofunction:: set_dihed

"""

import itertools
from types import MappingProxyType
from typing import List, Iterable, Union, Set, Optional, Type, Tuple, Mapping, Callable, Sequence
from collections import ChainMap

import numpy as np

import rdkit
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit.Chem import AllChem
from scm.plams.core.basejob import Job
from scm.plams import (Molecule, Atom, Bond, MoleculeError, add_to_class, Units,
                       Settings, AMSJob, ADFJob, Cp2kJob)

from .mol_split_cm import SplitMol
from .remove_atoms_cm import RemoveAtoms
from .optimize_rotmat import optimize_rotmat
from ..logger import logger
from ..workflows import WorkFlow, MOL, FORMULA, OPT
from ..mol_utils import fix_carboxyl, to_atnum
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file
from ..jobs import job_geometry_opt  # noqa: F401

__all__ = ['init_ligand_opt']

UFF = AllChem.UFFGetMoleculeForceField


def init_ligand_opt(ligand_df: SettingsDataFrame) -> None:
    """Initialize the ligand optimization procedure."""
    workflow = WorkFlow.from_template(ligand_df, name='ligand_opt')

    # Pull from the database; push unoptimized structures
    df_bool = workflow.from_db(ligand_df, FORMULA[0], read_mol=True)
    workflow.to_db(ligand_df, df_bool, columns=[MOL, FORMULA], status='no_opt')
    if not ligand_df.settings.optional.ligand.optimize:
        return None

    # Start the ligand optimization
    workflow(start_ligand_jobs, ligand_df, columns=[], index=df_bool[OPT])

    is_opt = (lig.properties.get('is_opt', False) for lig in ligand_df[MOL])
    ligand_df[OPT] = np.fromiter(is_opt, count=len(ligand_df), dtype=bool)

    # Push the optimized structures to the database
    workflow.to_db(ligand_df, df_bool, status='optimized')

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    mol_format = ligand_df.settings.optional.database.mol_format
    if mol_format:
        path = workflow.path
        mol_array = ligand_df.loc[df_bool[OPT], MOL].values
        mol_to_file(mol_array, path, mol_format=mol_format)


def _set_charge_adfjob(s: Settings, charge: int) -> None:
    s.input.charge = charge


def _set_charge_amsjob(s: Settings, charge: int) -> None:
    if s.get('uff') is not None:
        return
    s.input.ams.system.charge = charge


def _set_charge_cp2kjob(s: Settings, charge: int) -> None:
    if s.get('dft') is None:
        return
    s.input.force_eval.dft.charge = charge


ChargeFunc = Callable[[Settings, int], None]
CHARGE_FUNC_MAPPING: Mapping[Type[Job], ChargeFunc] = MappingProxyType({
    ADFJob: _set_charge_adfjob,
    AMSJob: _set_charge_amsjob,
    Cp2kJob: _set_charge_cp2kjob,
})


def start_ligand_jobs(ligand_list: Iterable[Molecule],
                      jobs: Iterable[Optional[Type[Job]]],
                      settings: Iterable[Optional[Settings]],
                      use_ff: bool = False, **kwargs) -> None:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform geometry optimizations."""
    _j1, job = jobs
    _s1, s = settings

    if job not in {None, AMSJob, ADFJob, Cp2kJob}:
        raise NotImplementedError(f"job1: {job.__class__} = {job!r}")
    elif use_ff:
        raise NotImplementedError(f"use_ff: {use_ff.__class__} = {use_ff!r}")

    if job is None:
        _start_ligand_jobs_uff(ligand_list)
    else:
        charge_func = CHARGE_FUNC_MAPPING[job]
        _start_ligand_jobs_plams(ligand_list, job, s, charge_func)
    return None


def _start_ligand_jobs_plams(ligand_list: Iterable[Molecule],
                             job: Type[Job], settings: Settings,
                             charge_func: ChargeFunc) -> None:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform geometry optimizations."""
    for ligand in ligand_list:
        try:
            optimize_ligand(ligand)
        except Exception as ex:
            logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)
            ligand.properties.is_opt = False
            continue

        s = Settings(settings)
        charge_func(s, int(sum(at.properties.get('charge', 0) for at in ligand)))
        ligand.job_geometry_opt(job, s, name='ligand_opt')
        ligand.round_coords()  # `.is_opt = True` is set by `ligand.job_geometry_opt()``
    return None


def _start_ligand_jobs_uff(ligand_list: Iterable[Molecule]) -> None:
    """Loop over all molecules in ``ligand_df.loc[idx]`` and perform geometry optimizations."""
    for ligand in ligand_list:
        logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization has started')
        try:
            optimize_ligand(ligand)
        except Exception as ex:
            ligand.properties.is_opt = False
            logger.error(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                         'has failed')
            logger.debug(f'{ex.__class__.__name__}: {ex}', exc_info=True)
        else:
            ligand.properties.is_opt = True
            logger.info(f'UFFGetMoleculeForceField: {ligand.properties.name} optimization '
                        'is successful')
    return None


def optimize_ligand(ligand: Molecule) -> None:
    """Optimize a ligand molecule."""
    anchor = ligand.properties.dummies

    # Split the branched ligand into linear fragments and optimize them individually
    bonds = split_mol(ligand, anchor)
    context = SplitMol(ligand, bonds)
    with context as mol_frags:
        cap_dict = ChainMap(*context._at_pairs)
        for mol in mol_frags:
            cap_list = [cap for at, cap in cap_dict.items() if at in mol]
            mol.set_dihed(180.0, anchor, cap_list)

    # Find the optimal dihedrals angle between the fragments
    for bond in bonds:
        modified_minimum_scan_rdkit(ligand, ligand.get_index(bond), anchor)

    # RDKit UFF can sometimes mess up the geometries of carboxylates: fix them
    fix_carboxyl(ligand)

    # Allign the ligand with the Cartesian X-axis.
    allign_axis(ligand, anchor)


def allign_axis(mol: Molecule, anchor: Atom):
    """Allign a molecule with the Cartesian X-axis; setting **anchor** as the origin."""
    try:
        idx = mol.atoms.index(anchor)
    except ValueError as ex:
        raise MoleculeError("The passed anchor is not in mol") from ex

    xyz = mol.as_array()  # Allign the molecule with the X-axis
    rotmat = optimize_rotmat(xyz, idx)
    xyz[:] = xyz@rotmat.T
    xyz -= xyz[idx]
    xyz[:] = xyz.round(decimals=3)
    mol.from_array(xyz)


def split_mol(plams_mol: Molecule, anchor: Atom) -> List[Bond]:
    """Split a molecule into multiple smaller fragments; returning the bonds that have to be broken.

    One fragment is created for every branch within **plams_mol**.

    Parameters
    ----------
    plams_mol : |plams.Molecule|
        The input molecule.

    anchor : |plams.Atom|
        An anchor atom which will be stored in the largest to-be returned fragment.

    Returns
    -------
    :class:`list` [|plams.Bond|]
        A list of plams bonds.

    """
    def _in_ring(bond: Bond) -> bool:
        """Check if one of the atoms in **bond** is part of a ring system."""
        return (plams_mol.in_ring(bond.atom1) or plams_mol.in_ring(bond.atom2))

    def _is_valid_bond(bond: Bond) -> bool:
        """Check if one atom in **bond** has at least 3 neighbours and the other at least 2."""
        n1, n2 = plams_mol.neighbors(bond.atom1), plams_mol.neighbors(bond.atom2)
        return (len(n1) >= 3 and len(n2) >= 2) or (len(n1) >= 2 and len(n2) >= 3)

    def _get_frag_size(bond: Bond) -> int:
        """Return the size of the largest fragment were **plams_mol** to be split along **bond**."""
        return plams_mol.get_frag_size(bond, anchor)

    # Temporary remove hydrogen atoms
    atom_gen = (at for at in plams_mol if at.atnum == 1)
    with RemoveAtoms(plams_mol, atom_gen):
        # Remove undesired bonds
        bond_list = [bond for bond in plams_mol.bonds if not _in_ring(bond)]

        # Remove even more undesired bonds
        for bond in reversed(bond_list):
            if not _is_valid_bond(bond):
                bond_list.remove(bond)

    atom_list = list(itertools.chain.from_iterable((bond.atom1, bond.atom2) for bond in bond_list))
    atom_set = {atom for atom in atom_list if atom_list.count(atom) >= 3}
    atom_dict = {atom: [bond for bond in atom.bonds if bond in bond_list] for atom in atom_set}

    # Fragment the molecule such that the anchor on the largest fragment
    ret = []
    for at, bond_list in atom_dict.items():
        # Can't directly iterate over bond_list as its size is modified
        iterator = range(len(bond_list[2:]))
        for _ in iterator:
            frag_size = [_get_frag_size(bond) for bond in bond_list]
            idx = np.argmax(frag_size)  # The index of the largest fragment
            bond = bond_list.pop(idx)
            ret.append(bond)
    return ret


@add_to_class(Molecule)
def get_frag_size(self, bond: Bond, anchor: Atom) -> int:
    """Return the size of the fragment containing **atom** if this instance was split into two by the breaking of **bond**.

    Parameters
    ----------
    bond : |plams.Bond|
        A PLAMS bond.

    anchor : |plams.Atom|
        A PLAMS atom. The size of the fragment containg this atom will be returned.

    Returns
    -------
    :class:`int`
        The number of atoms in the fragment containing **atom**.

    """  # noqa
    if bond not in self.bonds:
        raise MoleculeError('get_frag_size: The argument bond should be of type plams.Bond and '
                            'be part of the Molecule')
    elif anchor not in self.atoms:
        raise MoleculeError('get_frag_size: The argument atom should be of type plams.Atom and '
                            'be part of the Molecule')

    for at in self:
        at._visited = False

    frag1 = set()
    frag2 = set()

    def dfs(at1: Atom, atom_set: Set[Atom]):
        at1._visited = True
        atom_set.add(at1)
        for at2 in at1.neighbors():

            if at2._visited:
                continue
            dfs(at2, atom_set)

    bond.atom1._visited = bond.atom2._visited = True
    dfs(bond.atom1, frag1)
    dfs(bond.atom2, frag2)
    for at in self.atoms:
        del at._visited

    # fragment #1 contains **atom**
    if anchor in frag1:
        if bond.atom1 not in frag1:
            bond.atom1, bond.atom2 = bond.atom2, bond.atom1
        return len(frag1)

    if bond.atom1 not in frag2:
        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
    return len(frag2)


def get_dihed(atoms: Iterable[Atom], unit: str = 'degree') -> float:
    """Return the dihedral angle defined by a set of four atoms.

    Parameters
    ----------
    atoms : |Iterable| [|plams.atoms|]
        An iterable consisting of 4 PLAMS atoms

    unit : :class:`str`
        The output unit.

    Returns
    -------
    :class:`float`
        A dihedral angle expressed in **unit**.

    """
    at1, at2, at3, at4 = atoms
    vec1 = -np.array(at1.vector_to(at2))
    vec2 = np.array(at2.vector_to(at3))
    vec3 = np.array(at3.vector_to(at4))

    v1v2, v2v3 = np.cross(vec1, vec2), np.cross(vec3, vec2)
    v1v2_v2v3 = np.cross(v1v2, v2v3)
    v2_norm_v2 = vec2 / np.linalg.norm(vec2)
    epsilon = np.arctan2(v1v2_v2v3@v2_norm_v2, v1v2@v2v3)

    return Units.convert(epsilon, 'radian', unit)


@add_to_class(Molecule)
def neighbors_mod(self, atom: Atom, exclude: Union[int, str] = 1) -> List[Atom]:
    """A modified PLAMS function: Allows the exlucison of specific elements from the return list.

    Return a list of neighbors of **atom** within the molecule.
    Atoms with **atom** has to belong to the molecule.
    Returned list follows the same order as the **atom.bond** attribute.

    Parameters
    ----------
    atom : |plams.Atom|_
        The plams atom whose neighbours will be returned.

    exclude : |str|_ or |int|_
        Exclude all neighbours with a specific atomic number or symbol.

    Returns
    -------
    |list|_ [|plams.Atom|_]
        A list of all neighbours of **atom**.

    """
    exclude = to_atnum(exclude)
    if atom.mol != self:
        raise MoleculeError('neighbors: passed atom should belong to the molecule')
    return [b.other_end(atom) for b in atom.bonds if b.other_end(atom).atnum != exclude]


@add_to_class(Molecule)
def set_dihed(self, angle: float, anchor: Atom, cap: Sequence[Atom],
              opt: bool = True, unit: str = 'degree') -> None:
    """Change all valid dihedral angles into a specific value.

    Performs an inplace update of this instance.

    Parameters
    ----------
    angle : :class:`float`
        The desired dihedral angle.

    anchor : |plams.Atom|
        The ligand anchor atom.

    opt : :class:`bool`
        Whether or not the dihedral adjustment should be followed up by an RDKit UFF optimization.

    unit : :class:`str`
        The input unit.

    """
    cap_atnum = []
    for at in cap:
        cap_atnum.append(at.atnum)
        at.atnum = 0

    angle = Units.convert(angle, unit, 'degree')
    bond_iter = (bond for bond in self.bonds if bond.atom1.atnum != 1 and bond.atom2.atnum != 1
                 and bond.order == 1 and not self.in_ring(bond))

    # Correction factor for, most importantly, tri-valent anchors (e.g. P(R)(R)R)
    dihed_cor = angle / 2
    neighbors = anchor.neighbors()
    if len(neighbors) > 2:
        atom_list = [anchor] + sorted(neighbors, key=lambda at: -at.atnum)[:3]
        improper = get_dihed(atom_list)
        dihed_cor *= np.sign(improper)

    for bond in bond_iter:
        # Gather lists of all non-hydrogen neighbors
        n1, n2 = self.neighbors_mod(bond.atom1), self.neighbors_mod(bond.atom2)

        # Remove all atoms in `bond`
        n1 = [atom for atom in n1 if atom is not bond.atom2]
        n2 = [atom for atom in n2 if atom is not bond.atom1]

        # Remove all non-subsituted atoms
        # A special case consists of anchor atoms; they can stay
        if len(n1) > 1:
            n1 = [atom for atom in n1 if (len(self.neighbors_mod(atom)) > 1 or atom is anchor)]
        if len(n2) > 1:
            n2 = [atom for atom in n2 if (len(self.neighbors_mod(atom)) > 1 or atom is anchor)]

        # Set `bond` in an anti-periplanar conformation
        if n1 and n2:
            dihed = get_dihed((n1[0], bond.atom1, bond.atom2, n2[0]))
            if anchor not in bond:
                self.rotate_bond(bond, bond.atom1, angle - dihed, unit='degree')
            else:
                dihed -= dihed_cor
                self.rotate_bond(bond, bond.atom1, -dihed, unit='degree')
                dihed_cor *= -1

    for at, atnum in zip(cap, cap_atnum):
        at.atnum = atnum

    if opt:
        rdmol = molkit.to_rdmol(self)
        UFF(rdmol).Minimize()
        self.from_rdmol(rdmol)


def rdmol_as_array(rdmol: rdkit.Chem.Mol) -> np.ndarray:
    """Convert an rdkit molecule into an array of Cartesian coordinates."""
    def get_xyz(atom: rdkit.Chem.Atom) -> Tuple[float, float, float]:
        pos = conf.GetAtomPosition(atom.GetIdx())
        return (pos.x, pos.y, pos.z)

    conf = rdmol.GetConformer(id=-1)
    atoms = rdmol.GetAtoms()

    atom_count = len(atoms)
    count = 3 * atom_count
    shape = atom_count, 3

    iterator = itertools.chain.from_iterable(get_xyz(at) for at in atoms)
    ret = np.fromiter(iterator, count=count, dtype=float)
    ret.shape = shape
    return ret


def _find_idx(mol: Molecule, bond: Bond) -> List[int]:
    """Return the atomic indices of all atoms on the side of **bond.atom2**."""
    ret = []
    mol.set_atoms_id(start=0)
    for at in mol:
        at._visited = False

    def dfs(at1, mol):
        at1._visited = True
        ret.append(at1.id)
        for bond in at1.bonds:
            at2 = bond.other_end(at1)
            if not at2._visited:
                dfs(at2, mol)

    bond.atom1._visited = bond.atom2._visited = True
    dfs(bond.atom2, mol)

    mol.unset_atoms_id()
    return ret


def modified_minimum_scan_rdkit(ligand: Molecule, bond_tuple: Tuple[int, int],
                                anchor: Atom) -> None:
    """A modified version of the :func:`.global_minimum_scan_rdkit` function.

    * Uses the ligand vector as criteria rather than the energy.
    * Geometry optimizations are constrained during the conformation search.
    * Finish with a final unconstrained geometry optimization.

    See Also
    --------
    :func:`global_minimum_scan_rdkit<scm.plams.recipes.global_minimum.minimum_scan_rdkit>`:
        Optimize the molecule (RDKit UFF) with 3 different values for the given dihedral angle and
        find the lowest energy conformer.

        :param |Molecule| mol: The input molecule
        :param tuple bond_tuple: A 2-tuples containing the atomic indices of valid bonds
        :return |Molecule|: A copy of *mol* with a newly optimized geometry

    """
    # Define a number of variables and create 3 copies of the ligand
    angles = (-120, 0, 120)
    mol_list = [ligand.copy() for _ in range(3)]
    for angle, mol in zip(angles, mol_list):
        bond = mol[bond_tuple]
        atom = mol[bond_tuple[0]]
        mol.rotate_bond(bond, atom, angle, unit='degree')
    rdmol_list = [molkit.to_rdmol(mol, properties=False) for mol in mol_list]

    # Optimize the (constrained) geometry for all dihedral angles in angle_list
    # The geometry that yields the minimum energy is returned
    fixed = _find_idx(mol, bond)
    for rdmol in rdmol_list:
        ff = UFF(rdmol)
        for f in fixed:
            ff.AddFixedPoint(f)
        ff.Minimize()

    # Find the conformation with the optimal ligand vector
    cost_list = []
    try:
        i = mol.atoms.index(anchor)
    except ValueError:
        i = -1  # Default to the origin as anchor

    for rdmol in rdmol_list:
        xyz = rdmol_as_array(rdmol)
        if i == -1:  # Default to the origin as anchor
            xyz = np.vstack([xyz, [0, 0, 0]])
        rotmat = optimize_rotmat(xyz, i)
        xyz[:] = xyz@rotmat.T
        xyz -= xyz[i]
        cost = np.exp(xyz[:, 1:]).sum()
        cost_list.append(cost)

    # Perform an unconstrained optimization on the best geometry and update the geometry of ligand
    j = np.argmin(cost_list)
    rdmol_best = rdmol_list[j]
    UFF(rdmol).Minimize()
    ligand.from_rdmol(rdmol_best)
