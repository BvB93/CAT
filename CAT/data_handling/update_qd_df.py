"""A module for updating the index of quantum dot dataframes.

Index
-----
.. currentmodule:: CAT.data_handling.update_qd_df
.. autosummary::
    update_qd_df

API
---
.. autofunction:: update_qd_df

"""

import pandas as pd

from scm.plams import Molecule

from ..workflows import MOL, OPT, HDF5_INDEX

__all__ = ['update_qd_df']


def update_qd_df(qd_df: pd.DataFrame) -> None:
    """Update the index of **qd_df** with descriptors of the core(s) and ligand(s)."""
    # Construct the MultiIndex columns
    core = [_get_core_formula(mol) for mol in qd_df[MOL]]
    core_anchor = [None] * len(qd_df)
    ligand_smiles = [mol.properties.pop('ligand_smiles') for mol in qd_df[MOL]]
    ligand_anchor = [mol.properties.pop('ligand_anchor') for mol in qd_df[MOL]]

    # Reassign the index
    idx_iter = zip(*(core, core_anchor, ligand_smiles, ligand_anchor))
    names = ['core', 'core anchor', 'ligand smiles', 'ligand anchor']
    qd_df.index = pd.MultiIndex.from_tuples(idx_iter, names=names)

    # Set the ('opt', '') column to False
    qd_df[OPT] = False
    qd_df[HDF5_INDEX] = -1


def _get_core_formula(mol: Molecule) -> str:
    """Construct a string representing the molecular formula of the core."""
    # Create a dictionary with atomic symbols and matching atom counts
    core = {}
    for i in mol.properties.indices:
        at = mol[i]
        if 'anchor' not in at.properties:
            break  # The first ligand has been reached; abort and return
        symbol = at.symbol
        try:
            core[symbol] += 1
        except KeyError:
            core[symbol] = 1

    # Concatenate the dictionary into a single string
    return ''.join(f'{k}{v}' for k, v in sorted(core.items()))
