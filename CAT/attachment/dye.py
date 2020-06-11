"""A modules for combining (organic) molecules."""

from itertools import chain

import numpy as np
from scipy.spatial.distance import cdist

from rdkit.Chem import AllChem
from scm.plams import Molecule, Settings, to_rdmol

from .ligand_attach import rot_mol
from ..logger import logger

__all__ = ['label_lig', 'label_core', 'substitution', 'multi_substitution']

# flake8: noqa: N806


def ff_constrained_opt(mol, constrain=()):
    """Performs a constrained FF optimization on a PLAMS molecule.

    Optimisation with rdkit.Chem.AllChem.MMFFGetMoleculeForceField; a PLAMS molecule
    is converted to rdkit molecule, optimised with frozen atoms and converted back to
    PLAMS molecule.

    Parameters
    ----------
    mol : |plams.Molecule|
        A PLAMS molecule.

    constrain : list
        A list of indices of to-be frozen atoms.

    Returns
    -------
    |plams.Molecule|
        Optimized molecular structure

    """
    rdkit_mol = to_rdmol(mol)
    ff_type = AllChem.MMFFGetMoleculeForceField

    try:
        ff = ff_type(rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol))
        for f in constrain:
            ff.AddFixedPoint(f)
    except AttributeError:
        # It seems like the MMF forcefield is not available for all atoms (e.g. As)
        logger.warning(
            f"Failed to construct the {ff_type.__name__!r} forcefield for {mol.get_formula()}"
        )
        return mol

    ff.Minimize()
    mol.from_rdmol(rdkit_mol)

    return mol


def connect_ligands_to_core(lig_dict, core, user_min_dist):
    """Attaches multiple ligands to multiple copies of a single core.

    Returns a list of cores with attached ligands, each with the properties.min_distance attribute
    containing the smallest distance between ligand and core.

    Parameters
    ----------
    lig_dict : dict
        Ligands that are attached to the core molecule
        An iterable container consisting of a PLAMS molecules - ligand;
        A np.array of PLAMS atom coordinates - atom to be substituted;
        A vector array - bond between substituent and rest of the compound;
        An int - enumerating ligands
    core : |plams.Molecule|
        A PLAMS molecule - core on which ligands are attached
    user_min_dist : float
        Value for the minimal bond distance in the new molecule

    Returns
    -------
    list of PLAMS molecules
        A list of new molecule that are made core with different ligands attached

    """
    # Unpack the ligand dictionary
    lig_list = lig_dict['lig_list']
    lig_idx = lig_dict['lig_idx']
    lig_vec = lig_dict['lig_vec']
    lig_IDs = lig_dict['lig_ID']

    # Return if core.properties.vec has been exhausted
    try:
        core_vec = core.properties.vec[0]
    except IndexError:
        return []

    # Construct keyword arguments
    kwarg = get_args(core, lig_list, lig_idx)

    # Allign the ligands with the core; perform the rotation check
    lig_array, min_dist_array = rot_mol(lig_list, lig_vec, core_vec, **kwarg)

    # Reading list of ligands that are already attached to the core
    old_ligID = core.properties.get('ligID', '')

    # Combine the rotated ligands and core into new molecules
    ret = []
    for lig, xyz, min_dist, ligID in zip(lig_list, lig_array, min_dist_array, lig_IDs):
        # Copy and manipulate ligand
        lig_cp = lig.copy()
        lig_cp.guess_bonds()
        lig_cp.from_array(xyz)
        lig_cp.properties = Settings()

        # Adding new ligand to the list
        new_ligID = str(old_ligID) + str(ligID)

        # Copy and manipulate core
        core_h = core.properties.the_h
        core_cp = core.copy()
        core_cp.guess_bonds()
        core_cp.delete_atom(core_cp.closest_atom(core_h))
        core_cp.properties.ligID = new_ligID

        # Merge core and ligand
        lig_cp += core_cp
        lig_cp.properties.name = core.properties.name + "_" + lig.properties.name
        lig_cp.properties.min_distance = min_dist

        # Add bond between core and ligand
        the_h = lig_cp.closest_atom(lig_cp.properties.the_h)
        core_other = lig_cp.closest_atom(lig_cp.properties.coords_other_arrays[len(new_ligID)-1])
        lig_cp.add_bond(the_h, core_other)

        # Make list of atom indices to freeze, FrInd_rdkit, for optimization (rdkit starts at 0!)
        # Atom that ligand is conected to is not frozen
        Conn = list(lig_cp).index(core_other)
        FrInd_rdkit = [x for x in range(len(lig_cp)-len(core_cp), len(lig_cp)) if x != Conn]

        # FF for molecule with frozen core
        try:
            lig_cp = ff_constrained_opt(lig_cp, constrain=FrInd_rdkit)
        except ValueError:
            logger.warning("FF optimization error")

        # Distance between core and ligands
        FrInd_plams = (np.array(FrInd_rdkit)+1).tolist()     # plams counts from 1
        Fr_atoms = np.array([lig_cp[c].coords for c in FrInd_plams])
        OnlyLig = np.array(
                [lig_cp[i].coords for i in range(1, len(lig_cp)+1) if i not in FrInd_plams])
        min_dist = np.nanmin(cdist(OnlyLig, Fr_atoms))
        lig_cp.properties.min_distance = min_dist

        ret.append(lig_cp)

    return ret


def get_args(core, lig_list, lig_idx):
    """Extract the various arguments from core and ligand_list.

    Parameters
    ----------
    core : |plams.Molecule|
        A PLAMS molecule with plams_mol.properties attributes
    lig_list : list
        A list with PLAMS molecules with plams_mol.properties attributes
    lig_idx :
        A np.array of PLAMS atom coordinates - atom for substitution

    Returns
    -------
    dict
        A dictionary with properties needed for aligning ligand with the core

    """
    # Extract the various arguments from core and ligand_list
    core_other_atom = core.properties.coords_other_atom[0]
    lig_other = [lig[lig.properties.idx_other+1] for lig in lig_list]
    bond_length = np.array([core_other_atom.radius + lig.radius for lig in lig_other])

    # Roll through the arguments of core
    core.properties.the_h = core.properties.coords_h[0]
    at_h = core.closest_atom(core.properties.the_h)
    core.properties.vec = core.properties.vec[1:]
    core.properties.coords_other_atom = core.properties.coords_other_atom[1:]
    core.properties.coords_other = core.properties.coords_other[1:]
    core.properties.coords_h = core.properties.coords_h[1:]

    # Convert core into an array and mask core_h with np.nan
    core_array = core.as_array()
    idx = core.atoms.index(at_h)
    core_array[idx] = np.nan

    return {'core': core_array,
            'atoms_other': core_other_atom,
            'bond_length': bond_length,
            'dist_to_self': False,
            'idx': lig_idx,
            'ret_min_dist': True}


def label_core(mol):
    """Adds plams_mol.properties attribute to the core.

    Reads the atom indices from comment section in core's .xyz file and adds
    additional plams_mol.properties: coordinates of atom that will be substituted,
    bond vector between the substitution atom and its connection at the core,
    coordinates of the connection at the core

    Parameters
    ----------
    mol : |plams.Molecule|
        An input  PLAMS molecule with atom indices to be substituted in
        plams_mol.properties.comment

    Returns
    -------
    mol : |plams.Molecule|
        A PLAMS mol with additional plams_mol.properties

    """
    # Read the comment in the second line of the xyz file
    comment = mol.properties.comment
    comment = comment.split()
    mol.properties.core_len = len(mol)

    idx = np.array(comment, dtype=int)
    at_h = [mol[int(i)] for i in idx]
    at_other = [at.bonds[0].other_end(at) for at in at_h]

    dummy = Molecule()
    mol.properties.coords_h = dummy.as_array(atom_subset=at_h)
    mol.properties.coords_other = dummy.as_array(atom_subset=at_other)
    mol.properties.vec = mol.properties.coords_other - mol.properties.coords_h

    mol.properties.coords_h_atom = at_h
    mol.properties.coords_other_atom = at_other
    mol.properties.coords_other_arrays = [np.array(i.coords) for i in at_other]
    mol.guess_bonds()


def label_lig(mols):
    """Adds plams_mol.properties attribute to the ligand.

    Reads the atom index from comment section in ligand's .xyz file and adds additional
    plams_mol.properties: coordinates of atom that will be substituted; bond vector between
    the substitution atom and its connection at the ligand; coordinates of the connection
    at the ligand; ligands identity - ligands serial number

    Parameters
    ----------
    mol : |plams.Molecule|
        An input  PLAMS molecule with atom indices to be substituted in plams_mol.properties.comment

    Returns
    -------
    mol : |plams.Molecule|
        A PLAMS mol with additional plams_mol.properties

    """
    ligID = list(range(len(mols)))

    for mol, id in zip(mols, ligID):
        mol.properties.ligID = id
        mol.guess_bonds()

        # Read the comment in the second line of the xyz file
        comment = mol.properties.comment

        # Set an attirbute with indices and vectors
        mol.properties.vec = np.empty(3)

        # Fill the vector array, delete the marked atom in the ligand .xyz file
        at = mol[int(comment)]
        at_other = at.bonds[0].other_end(at)
        mol.properties.vec = np.array(at.coords) - np.array(at_other.coords)

        mol.delete_atom(at)
        mol.properties.idx_other = mol.atoms.index(at_other)
        mol.guess_bonds()


def substitution(input_ligands, input_cores, min_dist):
    """Substitutes atoms at the core with ligands.

    Parameters
    ----------
    input_ligands : list
        A list of input PLAMS molecules with index of an atom to be substituted in
            plams_mol.properties.comment input_cores : list
        A list of input PLAMS molecules with list of indices of atoms to be
            substituted in plams_mol.properties.comment min_dist : float
        Minimal distance between core and attached ligands

    Returns
    -------
    list
        New molecules that are made of ligands attached to the core at position
        of the first index in the list
    """
    lig_idx = np.array([lig.properties.idx_other for lig in input_ligands])
    lig_vec = np.array([lig.properties.vec for lig in input_ligands])
    lig_ID = [lig.properties.ligID for lig in input_ligands]
    lig_dict = {'lig_list': input_ligands, 'lig_idx': lig_idx, 'lig_vec': lig_vec,
                'lig_ID': lig_ID}
    ret = (connect_ligands_to_core(lig_dict, core, min_dist) for core in input_cores)
    return list(chain.from_iterable(ret))


def multi_substitution(input_ligands, input_cores, n=1):
    """Attach ligands to cores; repeat this process n times."""
    ret = []
    mol_list = input_cores
    for _ in range(n):
        mol_list = substitution(input_ligands, mol_list)
        ret.append(mol_list)
    return list(chain.from_iterable(ret))
