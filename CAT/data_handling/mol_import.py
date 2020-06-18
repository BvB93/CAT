"""A module related to the importing of molecules.

Index
-----
.. currentmodule:: CAT.data_handling.mol_import
.. autosummary::
    read_mol
    read_mol_xyz
    read_mol_pdb
    read_mol_mol
    read_mol_smiles
    read_mol_plams
    read_mol_rdkit
    read_mol_folder
    read_mol_txt
    get_charge_dict
    set_mol_prop
    set_atom_prop
    set_qd
    print_exception

API
---
.. autofunction:: read_mol
.. autofunction:: read_mol_xyz
.. autofunction:: read_mol_pdb
.. autofunction:: read_mol_mol
.. autofunction:: read_mol_smiles
.. autofunction:: read_mol_plams
.. autofunction:: read_mol_rdkit
.. autofunction:: read_mol_folder
.. autofunction:: read_mol_txt
.. autofunction:: get_charge_dict
.. autofunction:: set_mol_prop
.. autofunction:: set_atom_prop
.. autofunction:: set_qd
.. autofunction:: print_exception

"""

import os
import itertools
from types import MappingProxyType
from string import ascii_letters
from typing import Dict, Iterable, List, Sequence, Optional, Mapping, Callable

import numpy as np

from scm.plams import (Molecule, Atom, Settings)
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem, RDLogger

from ..logger import logger
from ..data_handling.validate_mol import validate_mol

__all__ = ['read_mol', 'set_mol_prop']

# Supress the rdkit logger
_logger = RDLogger.logger()
_logger.setLevel(RDLogger.CRITICAL)


def read_mol(input_mol: Iterable[Settings]) -> List[Molecule]:
    """Checks the filetypes of the input molecules.

    Sets the molecules' properties and returns a list of plams molecules.

    Parameters
    ----------
    input_mol : |list|_ [|Settings|_]
        An iterable consisting of dictionaries with input settings per mol.

    Returns
    -------
    |plams.Molecule|_
        A list of plams Molecules.

    """
    # Create a list of PLAMS molecules
    mol_list = []
    append = mol_list.append
    extend = mol_list.extend

    for mol_dict in input_mol:
        try:
            read_mol = EXTENSION_MAPPING[mol_dict.type]
        except KeyError as ex:
            logger.error(f'{ex.__class__.__name__}: {ex}')
            continue

        mol = read_mol(mol_dict)
        if not mol:  # Failed to import any molecules
            continue

        if isinstance(mol, list):  # if mol is a list of molecules
            extend(mol)

        elif mol_dict.is_qd:  # A quantum dot molecule
            set_qd(mol, mol_dict)
            append(mol)

        else:  # A core or ligand molecule
            if mol_dict.guess_bonds:
                mol.guess_bonds()
            set_mol_prop(mol, mol_dict)
            append(mol)

    return mol_list


def read_mol_xyz(mol_dict: Settings) -> Optional[Molecule]:
    """Read an .xyz file."""
    try:
        mol = Molecule(mol_dict.mol, inputformat='xyz')
        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        if not mol_dict.is_core:
            canonicalize_mol(mol)
        return mol
    except Exception as ex:
        print_exception('read_mol_xyz', mol_dict.name, ex)


def read_mol_pdb(mol_dict: Settings) -> Optional[Molecule]:
    """Read a .pdb file."""
    try:
        mol = molkit.readpdb(mol_dict.mol)
        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        if not mol_dict.is_core:
            canonicalize_mol(mol)
        return mol
    except Exception as ex:
        print_exception('read_mol_pdb', mol_dict.name, ex)


def read_mol_mol(mol_dict: Settings) -> Optional[Molecule]:
    """Read a .mol file."""
    try:
        mol = molkit.from_rdmol(Chem.MolFromMolFile(mol_dict.mol, removeHs=False))
        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        if not mol_dict.is_core:
            canonicalize_mol(mol)
        return mol
    except Exception as ex:
        print_exception('read_mol_mol', mol_dict.name, ex)


def read_mol_smiles(mol_dict: Settings) -> Optional[Molecule]:
    """Read a SMILES string."""
    try:
        mol = _from_smiles(mol_dict.mol)
        mol.properties.smiles = Chem.CanonSmiles(mol_dict.mol)

        if mol_dict.get('indices'):
            for i in mol_dict.indices:
                mol[i].properties.anchor = True

        canonicalize_mol(mol)
        if mol_dict.get('indices'):
            mol_dict.indices = tuple(i for i, at in enumerate(mol, 1) if
                                     at.properties.pop('anchor', False))

        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        return mol
    except Exception as ex:
        print_exception('read_mol_smiles', mol_dict.name, ex)


def _from_smiles(smiles: str) -> Molecule:
    sanitize = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    _rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.rdmolops.SanitizeMol(_rdmol, sanitizeOps=sanitize)

    rdmol = Chem.AddHs(_rdmol)
    rdmol.SetProp('smiles', smiles)
    return molkit.get_conformations(rdmol, 1, None, None, rms=0.1)


def read_mol_plams(mol_dict: Settings) -> Optional[Molecule]:
    """Read a PLAMS molecule."""
    try:
        mol = mol_dict.mol
        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        if not mol_dict.is_core:
            canonicalize_mol(mol)
        return mol
    except Exception as ex:
        print_exception('read_mol_plams', mol_dict.name, ex)


def read_mol_rdkit(mol_dict: Settings) -> Optional[Molecule]:
    """Read a RDKit molecule."""
    try:
        mol = molkit.from_rdmol(mol_dict.mol)
        if mol_dict.guess_bonds and not mol_dict.is_qd:
            mol.guess_bonds()
        if not mol_dict.is_core:
            canonicalize_mol(mol)
        return mol
    except Exception as ex:
        print_exception('read_mol_rdkit', mol_dict.name, ex)


def read_mol_folder(mol_dict: Settings) -> Optional[List[Molecule]]:
    """Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder."""
    try:
        mol_type = 'input_cores' if mol_dict.is_core else 'input_ligands'

        _file_list = os.listdir(mol_dict.mol)
        optional_dict = Settings({k: v for k, v in mol_dict.items() if k not in ('mol', 'path')})
        file_list = [{i: optional_dict} for i in _file_list]

        validate_mol(file_list, mol_type, mol_dict.path)
        return read_mol(file_list)
    except Exception as ex:
        print_exception('read_mol_folder', mol_dict.name, ex)


def read_mol_txt(mol_dict: Settings) -> Optional[List[Molecule]]:
    """Read a plain text file containing one or more SMILES strings."""
    try:
        row = 0 if 'row' not in mol_dict else mol_dict.row
        column = 0 if 'column' not in mol_dict else mol_dict.column
        mol_type = 'input_cores' if mol_dict.is_core else 'input_ligands'

        with open(mol_dict.mol, 'r') as f:
            iterator = itertools.islice(f, row, None)
            _file_list = [i.rstrip('\n').split()[column] for i in iterator if i]
        optional_dict = Settings({k: v for k, v in mol_dict.items() if k not in ('mol', 'path')})
        file_list = [{i: optional_dict} for i in _file_list]

        validate_mol(file_list, mol_type, mol_dict.path)
        return read_mol(file_list)
    except Exception as ex:
        print_exception('read_mol_txt', mol_dict.name, ex)


EXTENSION_MAPPING: Mapping[str, Callable[[Settings], Optional[Molecule]]] = MappingProxyType({
    'xyz': read_mol_xyz,
    'pdb': read_mol_pdb,
    'mol': read_mol_mol,
    'smiles': read_mol_smiles,
    'folder': read_mol_folder,
    'txt': read_mol_txt,
    'plams_mol': read_mol_plams,
    'rdmol': read_mol_rdkit
})


def canonicalize_mol(mol: Molecule, inplace: bool = True) -> Optional[Molecule]:
    """Take a PLAMS molecule and sort its atoms based on their canonical rank.

    .. _rdkit.Chem.CanonicalRankAtoms: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.CanonicalRankAtoms

    Examples
    --------
    .. code:: python

        >>> from scm.plams import Molecule, from_smiles

        # Methane
        >>> mol: Molecule = from_smiles('C')
        >>> print(mol)  # doctest: +SKIP
        Atoms:
            1         H      0.640510      0.640510     -0.640510
            2         H      0.640510     -0.640510      0.640510
            3         C      0.000000      0.000000      0.000000
            4         H     -0.640510      0.640510      0.640510
            5         H     -0.640510     -0.640510     -0.640510

        >>> canonicalize_mol(mol)
        >>> print(mol)  # doctest: +SKIP
        Atoms:
            1         C      0.000000      0.000000      0.000000
            2         H     -0.640510     -0.640510     -0.640510
            3         H     -0.640510      0.640510      0.640510
            4         H      0.640510     -0.640510      0.640510
            5         H      0.640510      0.640510     -0.640510

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    inplace : bool
        If ``True``, perform an inplace update of **mol** rather than returning
        a new :class:`Molecule` instance.

    Returns
    -------
    |plams.Molecule|_
        Optional: if ``inplace=False``, return a copy of **mol** with its atoms sorted by their
        canonical rank.

    See Also
    --------
    * rdkit.Chem.CanonicalRankAtoms_: Returns the canonical atom ranking for each atom of a
      molecule fragment.

    """  # noqa
    rdmol = molkit.to_rdmol(mol)
    idx_collection = Chem.CanonicalRankAtoms(rdmol)

    # Reverse sort Molecule.atoms by the atomic indices in idx_collection
    if inplace:
        mol.atoms = [at for _, at in sorted(zip(idx_collection, mol.atoms), reverse=True)]
        return
    else:
        ret = mol.copy()
        ret.atoms = [at for _, at in sorted(zip(idx_collection, ret.atoms), reverse=True)]
        return ret


def get_charge_dict() -> Dict[str, int]:
    """Create a dictionary of elements and their formal atomic charge."""
    # Create a list of atomic charges and elements
    charges = (1, 2, -3, -2, -1, 2)
    elements = (
        ('Li', 'Na', 'K', 'Rb', 'Cs'),      # Group 1: Alkali metals
        ('Be', 'Mg', 'Ca', 'Sr', 'Ba'),     # Group 2: Alkaline earth metals
        ('N', 'P', 'As', 'Sb', 'Bi'),       # Group 15: pnictogens
        ('O', 'S', 'Se', 'Te', 'Po'),       # Group 16: Chalcogens
        ('H', 'F', 'Cl', 'Br', 'I', 'At'),  # Group 17: Halogens
        ('Cd', 'Pb')                        # Misc
    )

    # Combine the elements and atomic charges into a dictionary
    values = itertools.chain.from_iterable([i] * len(j) for i, j in zip(charges, elements))
    keys = itertools.chain(*elements)

    return dict(zip(keys, values))


#: A dictionary with default atomic charges
CHARGE_MAPPING: Mapping[str, int] = MappingProxyType(get_charge_dict())


def set_mol_prop(mol: Molecule, mol_dict: Settings) -> None:
    """Set molecular and atomic properties."""
    if mol_dict.is_core:
        residue_name = 'COR'
        mol.properties.name = mol.get_formula()
    else:
        residue_name = 'LIG'
        mol.properties.name = mol_dict.name

    mol.properties.dummies = mol_dict.indices
    mol.properties.path = mol_dict.path
    mol.properties.job_path = []

    # Prepare a generator of letters for pdb_info.Name
    alphabet = itertools.combinations(ascii_letters, 2)

    # Set the atomic properties
    for atom, i in zip(mol, itertools.cycle(alphabet)):
        set_atom_prop(atom, i, residue_name)

    if not mol.properties.smiles:
        mol.properties.smiles = Chem.MolToSmiles(
            Chem.RemoveHs(molkit.to_rdmol(mol)), canonical=True
        )


def set_atom_prop(atom: Atom, at_id: Sequence[str], residue_name: str) -> None:
    """Set atomic properties."""
    symbol = '{:4}'.format(atom.symbol + ''.join(at_id))

    # Add a number of properties to atom
    atom.properties.pdb_info.ResidueName = residue_name
    atom.properties.pdb_info.Occupancy = 1.0
    atom.properties.pdb_info.TempFactor = 0.0
    atom.properties.pdb_info.ResidueNumber = 1
    atom.properties.pdb_info.Name = symbol
    atom.properties.pdb_info.ChainId = 'A'

    # Changes hydrogen and carbon from heteroatom to atom
    if atom.symbol in ('H', 'C'):
        atom.properties.pdb_info.IsHeteroAtom = False
    else:
        atom.properties.pdb_info.IsHeteroAtom = True

    # Sets the formal atomic charge
    if atom.properties.charge:
        return

    # Default to a charge of 0 if no charge is available for that specific element
    if atom.symbol not in CHARGE_MAPPING:
        atom.properties.charge = 0
        return

    # Update the charge of non-hypervalent atoms
    total_bonds = int(sum([bond.order for bond in atom.bonds]))
    default_charge = CHARGE_MAPPING[atom.symbol]
    abs_charge = abs(default_charge)
    sign = -1 * int(default_charge / abs_charge)

    # Take the default charge and correct for the number (and order) of bonds
    atom.properties.charge = default_charge + sign * total_bonds
    if total_bonds <= abs_charge:
        return

    # Update formal atomic charges for hypervalent atoms
    if total_bonds is 2 + abs_charge:
        atom.properties.charge -= 2 * sign
    elif total_bonds is 4 + abs_charge:
        atom.properties.charge -= 4 * sign
    elif total_bonds >= 6 + abs_charge:
        atom.properties.charge -= 6 * sign
    return


def set_qd(qd: Molecule, mol_dict: Settings) -> Molecule:
    """Update quantum dots imported by :func:`.read_mol`."""
    # Create ligand (and anchor) molecules
    ligand = molkit.from_smiles(mol_dict.ligand_smiles)
    ligand_rdmol = molkit.to_rdmol(ligand)
    anchor = molkit.from_smiles(mol_dict.ligand_anchor)
    anchor_rdmol = molkit.to_rdmol(anchor)
    qd_rdmol = molkit.to_rdmol(qd)

    # Create arrays of atomic indices of the core and ligands
    lig_idx = 1 + np.array(qd_rdmol.GetSubstructMatches(ligand_rdmol))
    core_idx = np.arange(1, len(qd))[~lig_idx]
    lig_idx = lig_idx.ravel().tolist()
    core_idx = core_idx.tolist()

    # Guess bonds
    if mol_dict.guess_bonds:
        qd.guess_bonds(atom_subset=[qd[i] for i in lig_idx])

    # Reorder all atoms: core atoms first followed by ligands
    qd.atoms = [qd[i] for i in core_idx] + [qd[j] for i in lig_idx for j in i]

    # Construct a list with the indices of all ligand anchor atoms
    core_idx_max = 1 + len(core_idx)
    _anchor_idx = ligand_rdmol.GetSubstructMatch(anchor_rdmol)[0]
    start = core_idx_max + _anchor_idx
    stop = core_idx_max + _anchor_idx + np.product(lig_idx.shape)
    step = len(ligand)
    anchor_idx = list(range(start, stop, step))

    # Update the properties of **qd**
    for i in anchor_idx:
        qd[i].properties.anchor = True
    qd.properties.indices = list(range(1, core_idx_max)) + anchor_idx
    qd.properties.job_path = []
    qd.properties.name = mol_dict.name
    qd.properties.path = mol_dict.path
    qd.properties.ligand_smiles = Chem.CanonSmiles(mol_dict.ligand_smiles)
    qd.properties.ligand_anchor = f'{ligand[_anchor_idx].symbol}{_anchor_idx}'

    # Update the pdb_info of all atoms
    for i, at in enumerate(qd, 1):
        at.properties.pdb_info.SerialNumber = i
        if i <= core_idx_max:  # A core atom
            at.properties.pdb_info.ResidueNumber = 1
        else:  # A ligand atom
            at.properties.pdb_info.ResidueNumber = 2 + int((i - core_idx_max) / len(ligand))


def print_exception(func_name: str, mol_name: str, ex: Exception) -> None:
    """Manages the printing of exceptions upon failing to import a molecule."""
    ex_name = ex.__class__.__name__
    err = f'CAT.{func_name}() failed to parse {repr(mol_name)}; Caught exception: {ex_name}'
    logger.error(err)
