"""A module for attaching multiple non-unique ligands to a single quantum dot."""

import os
from typing import (Iterable, Any, overload, Sequence, MutableSequence,
                    Collection, List, Union)

import numpy as np
import pandas as pd

from rdkit import Chem
from scm.plams import Molecule, MoleculeError

from .workflows import WorkFlow
from .mol_utils import to_symbol
from .data_handling import mol_to_file
from .data_handling.mol_import import read_mol
from .data_handling.validate_mol import validate_mol
from .attachment.ligand_opt import optimize_ligand, allign_axis
from .attachment.ligand_anchoring import find_substructure
from .attachment.ligand_attach import ligand_to_qd

__all__ = ['init_multi_ligand']


def init_multi_ligand(qd_df):
    """Initialize the multi-ligand attachment procedure."""
    workflow = WorkFlow.from_template(qd_df, name='multi_ligand')

    if workflow.dummy is not None:
        sequence = [to_symbol(i) for i in workflow.dummy]
    elif workflow.f is not None:
        sequence = [str(i) for i in workflow.f]
    else:
        raise TypeError("'workflow.f' and 'workflow.dummy' cannot be both 'None'")

    columns_iter1 = ('/'.join(item for item in sequence[:i]) for i in range(1, 1+len(sequence)))
    columns_iter2 = (('multi ligand', i) for i in columns_iter1)
    columns = pd.MultiIndex.from_tuples(columns_iter2, names=qd_df.columns.names)
    for k in columns:
        qd_df[k] = None

    workflow(multi_lig, qd_df, columns=columns, no_loc=True)

    if workflow.mol_format is not None:
        path = workflow.path
        mol_format = workflow.mol_format
        mol_ar_flat = qd_df[columns].values.ravel()
        mol_to_file(mol_ar_flat, path, mol_format=mol_format)


@overload
def multi_lig(qd_series: pd.Series, ligands: Iterable[str],
              dummy: Sequence[Union[str, int]], f: None,
              **kwargs: Any) -> pd.DataFrame:
    ...
@overload  # noqa: E302
def multi_lig(qd_series: pd.Series, ligands: Iterable[str],
              dummy: None, f: Sequence[float],
              **kwargs: Any) -> pd.DataFrame:
    ...
def multi_lig(qd_series, ligands, dummy=None, f=None, **kwargs):  # noqa: E302
    """Attach multiple non-unique **ligands** to each qd in **qd_series**."""
    # Read and parse the SMILES strings
    ligands = smiles_to_lig(list(ligands),
                            functional_groups=kwargs['functional_groups'],
                            opt=kwargs['opt'],
                            split=kwargs['split'])

    mol_format = kwargs.get('mol_format')
    _path = kwargs.get('path')
    if mol_format is not _path is not None:
        path = os.path.join(str(_path).rsplit(os.sep, 1)[0], 'ligand')
        mol_to_file(ligands, path, mol_format=mol_format)

    if f is not None:
        raise NotImplementedError("'f != None' is not yet implemented")

    if dummy is not None:
        return _multi_lig_dummy(qd_series, ligands, kwargs['path'], dummy, kwargs['allignment'])
    elif f is not None:
        return [[NotImplemented]]
    else:
        raise TypeError("'f' and 'dummy' cannot be both 'None'")


def _multi_lig_dummy(qd_series, ligands, path, dummy, allignment) -> np.ndarray:
    """Gogogo."""
    ret = np.empty((len(ligands), len(qd_series)), dtype=object)
    for i, qd in enumerate(qd_series):
        qd = qd.copy()

        for j, (ligand, atnum) in enumerate(zip(ligands, dummy)):
            try:
                atoms = [at for at in qd if at.atnum == atnum]
                assert atoms
            except AssertionError as ex:
                raise MoleculeError(f'Failed to identify {to_symbol(atnum)!r} in '
                                    f'{qd.get_formula()!r}') from ex
            else:
                for at in atoms:
                    qd.delete_atom(at)

            coords = Molecule.as_array(None, atom_subset=atoms)
            qd.properties.dummies = np.array(coords, ndmin=2, dtype=float)
            qd = ligand_to_qd(qd, ligand, path=path,
                              allignment=allignment,
                              idx_subset=qd.properties.indices)
            ret[j, i] = qd
    return ret


def _multi_lig_f(qd_series, ligands, path, f, **kwargs):
    raise NotImplementedError


def smiles_to_lig(smiles: MutableSequence[str],
                  functional_groups: Collection[Chem.Mol],
                  opt: bool = True, split: bool = True) -> List[Molecule]:
    """Parse and convert all **smiles** strings into Molecules."""
    # Convert the SMILES strings into ligands
    validate_mol(smiles, 'input_ligands')
    ligands = [find_substructure(lig, functional_groups, split)[0] for lig in read_mol(smiles)]

    # Optimize the ligands
    process_mol = optimize_ligand if opt else lambda n: allign_axis(n, n.properties.dummies)
    for lig in ligands:
        process_mol(lig)
    return ligands
