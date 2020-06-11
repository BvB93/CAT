"""A module designed for finding ligand functional groups.

Index
-----
.. currentmodule:: CAT.attachment.ligand_anchoring
.. autosummary::
    init_ligand_anchoring
    get_functional_groups
    _smiles_to_rdmol
    find_substructure
    substructure_split
    _get_df

API
---
.. autofunction:: init_ligand_anchoring
.. autofunction:: get_functional_groups
.. autofunction:: find_substructure
.. autofunction:: _smiles_to_rdmol
.. autofunction:: substructure_split
.. autofunction:: _get_df

"""

from itertools import chain
from typing import Sequence, List, Tuple, Optional, Iterable, Callable

import pandas as pd

from scm.plams import Molecule, Settings
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from ..logger import logger
from ..utils import get_template
from ..mol_utils import separate_mod   # noqa: F401
from ..workflows import MOL, FORMULA, HDF5_INDEX, OPT
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.validate_mol import santize_smiles

__all__ = ['init_ligand_anchoring']


def init_ligand_anchoring(ligand_df: SettingsDataFrame) -> SettingsDataFrame:
    """Initialize the ligand functional group searcher.

    Parameters
    ----------
    ligand_df : |CAT.SettingsDataFrame|_
        A dataframe of valid ligands.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Unpack arguments
    settings = ligand_df.settings.optional
    split = settings.ligand.split
    functional_groups = settings.ligand.functional_groups

    # Find all functional groups; return a copy of each mol for each functional group
    mol_list = []
    for lig in ligand_df[MOL]:
        # Functional group search
        dummies = lig.properties.dummies
        if not dummies:
            mol_list += find_substructure(lig, functional_groups, split)
            continue

        # Manual specification of a functional group
        if len(dummies) == 1:  # optional.ligand.split = False
            lig.properties.dummies = dummies[0] - 1
            _split = False
        elif len(dummies) == 2:  # optional.ligand.split = True
            lig.properties.dummies = tuple(i - 1 for i in dummies)
            _split = True

        mol_list += [substructure_split(lig, lig.properties.dummies, split=_split)]

    # Convert the results into a dataframe
    return _get_df(mol_list, ligand_df.settings)


def _get_df(mol_list: Sequence[Molecule],
            settings: Settings) -> SettingsDataFrame:
    """Create and return a new ligand dataframe.

    Parameters
    ----------
    mol_list : |list|_ [|plams.Molecule|_]
        A list of PLAMS molecules.

    settings : |Settings|_
        A Settings instance containing all CAT parameters.

    Returns
    -------
    |CAT.SettingsDataFrame|_
        A dataframe of ligands with functional groups that can serve as valid anchor points.

    """
    # Create the dataframe index and columns
    idx_tuples = [(mol.properties.smiles, mol.properties.anchor) for mol in mol_list]
    idx = pd.MultiIndex.from_tuples(idx_tuples, names=['smiles', 'anchor'])
    columns_tuples = [MOL, FORMULA, HDF5_INDEX, OPT]
    columns = pd.MultiIndex.from_tuples(columns_tuples, names=['index', 'sub index'])

    # Create, fill and return the dataframe
    df = SettingsDataFrame(-1, index=idx, columns=columns, settings=settings)
    df[MOL] = mol_list
    df[FORMULA] = [lig.get_formula() for lig in df[MOL]]
    df[OPT] = False
    return df[~df.index.duplicated(keep='first')]  # Remove duplicate indices


def get_functional_groups(functional_groups: Optional[Iterable[str]] = None,
                          split: bool = True) -> Tuple[Chem.Mol, ...]:
    """Construct a list of RDKit molecules representing functional groups.

    Parameters
    ----------
    functional_groups : |list|_ [|str|_]
        Optional: A list of smiles strings representing functional groups.
        Will default to templates provided by CAT if ``None``.

    split : bool
        If templates should be pulled from the ``['split']`` or ``['no_split']`` block.
        Only relevant if **functional_groups** is ``None``.

    Returns
    -------
    |tuple|_ [|Chem.Mol|_]
        A list of RDKit molecules constructed from either **functional_group** or
        the default smiles templates in CAT.

    """
    # The user has, explicitly, provided functional groups
    if functional_groups:
        return tuple(_smiles_to_rdmol(smiles) for smiles in functional_groups)

    # Read functional groups from the default CAT SMILES templates
    if split:
        func_groups = get_template('smiles.yaml').split
    else:
        func_groups = get_template('smiles.yaml').no_split

    return tuple(_smiles_to_rdmol(smiles) for smiles in func_groups)


def _smiles_to_rdmol(smiles: str) -> Chem.Mol:
    """Convert a SMILES string into an rdkit Mol; supports explicit hydrogens."""
    # RDKit tends to remove explicit hydrogens if SANITIZE_ADJUSTHS is enabled
    sanitize = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.rdmolops.SanitizeMol(mol, sanitizeOps=sanitize)
    except Exception as ex:
        raise ex.__class__(f'Failed to parse the following SMILES string: {repr(smiles)}\n\n{ex}')
    return mol


def find_substructure(ligand: Molecule,
                      func_groups: Iterable[Chem.Mol],
                      split: bool = True,
                      condition: Optional[Callable[[int], bool]] = None) -> List[Molecule]:
    """Identify interesting functional groups within the ligand.

    Parameters
    ----------
    ligand : |plams.Molecule|_
        The ligand molecule.

    func_groups : |tuple|_ [|Chem.Mol|_]
        A collection of RDKit molecules representing functional groups.

    split : bool
        If a functional group should be split from **ligand** (``True``) or not (``False``).

    Returns
    -------
    |list|_ [|plams.Molecule|_]
        A list of ligands.
        A single copy of **ligand** is created for each identified functional group,
        removing parts of the functional group if required (see **split**).
        An empty list is returned if no valid functional groups are found.

    """
    rdmol = molkit.to_rdmol(ligand)

    # Searches for functional groups (defined by functional_group_list) within the ligand
    get_match = rdmol.GetSubstructMatches
    matches = chain.from_iterable(get_match(mol, useChirality=True) for mol in func_groups)

    # Remove all duplicate matches, each heteroatom (match[0]) should have <= 1 entry
    ligand_indices = []
    ref = []
    for idx_tup in matches:
        i, *_ = idx_tup
        if i in ref:
            continue  # Skip duplicates

        ligand_indices.append(idx_tup)
        ref.append(i)

    if condition is not None:
        if not condition(len(ligand_indices)):
            err = (f"Failed to satisfy the passed condition ({condition!r}) for "
                   f"ligand: {ligand.properties.name!r}")
            logger.error(err)
            return []
    if ligand_indices:
        return [substructure_split(ligand, tup, split) for tup in ligand_indices]
    else:
        err = (f"No functional groups were found (optional.ligand.split = {split!r}) for "
               f"ligand: {ligand.properties.name!r}")
        logger.error(err)
        return []


def substructure_split(ligand: Molecule,
                       idx: Tuple[int, int],
                       split: bool = True) -> Molecule:
    """Delete the hydrogen or mono-/polyatomic counterion attached to the functional group.

    Sets the charge of the remaining heteroatom to -1 if ``split=True``.

    Parameters
    ----------
    ligand: |plams.Molecule|_
        The ligand molecule.

    idx : |tuple|_ [|int|_]
        A tuple with 2 atomic indices associated with a functional group.

    split : bool
        If a functional group should be split from **ligand** (``True``) or not (``False``).

    Returns
    -------
    |plams.Molecule|_
        A copy of **ligand**, with part of its functional group removed (see **split**).

    """
    lig = ligand.copy()
    at1 = lig[idx[0] + 1]
    at2 = lig[idx[-1] + 1]

    if split:
        lig.delete_atom(at2)
        mol_list = lig.separate_mod()
        for mol in mol_list:
            if at1 not in mol:
                continue

            lig = mol
            break

        # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
        if not at1.properties.charge or at1.properties.charge == 0:
            at1.properties.charge = -1

    # Update ligand properties
    lig.properties.dummies = at1
    lig.properties.anchor = at1.symbol + str(lig.atoms.index(at1) + 1)
    lig.properties.charge = sum(atom.properties.charge for atom in lig if atom.properties.charge)

    # Update the ligand smiles string
    rdmol = molkit.to_rdmol(lig)
    smiles = Chem.MolToSmiles(rdmol)
    lig.properties.smiles = Chem.CanonSmiles(smiles)
    lig.properties.name = santize_smiles(lig.properties.smiles) + '@' + lig.properties.anchor
    lig.properties.path = ligand.properties.path

    return lig
