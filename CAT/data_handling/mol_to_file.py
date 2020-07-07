"""A module designed for exporting PLAMS :class:`Molecule` instances to common file formats.

Index
-----
.. currentmodule:: CAT.data_handling.mol_to_file
.. autosummary::
    mol_to_file

API
---
.. autofunction:: mol_to_file

"""

import os
from types import MappingProxyType
from typing import Iterable, Optional, Container, Mapping, Union, Callable
from os.path import join, isdir, isfile, exists

from scm.plams import Molecule, writepdb

__all__ = ['mol_to_file']

MolExportFunc = Callable[[Molecule, Union[str, bytes, os.PathLike]], None]

#: A mapping of file extensions to a Callable for Molecule exporting
EXPORT_MAPPING: Mapping[str, MolExportFunc] = MappingProxyType({
    'pdb': writepdb,
    'xyz': Molecule.write,
    'mol': Molecule.write,
    'mol2': Molecule.write
})


def mol_to_file(mol_list: Iterable[Molecule],
                path: Union[None, str, 'os.PathLike[str]'] = None,
                overwrite: bool = True,
                mol_format: Container[str] = frozenset({'xyz', 'pdb'})) -> None:
    """Export all molecules in **mol_list** to .pdb, .xyz, .mol or .mol2 files.

    Parameters
    ----------
    mol_list: |list|_ [|plams.Molecule|_]
        An iterable consisting of PLAMS molecules.

    path : str
        Optional: The path to the directory where the molecules will be stored.
        Defaults to the current working directory if ``None``.

    overwrite : bool
        If previously generated files can be overwritten or not.

    mol_format : |list|_ [|str|_]
        A list of strings with the to-be exported file types.
        Accepted values are ``"xyz"``, ``"pdb"``, ``"mol"`` and/or ``"mol2"``.

    Raises
    ------
    FileNotFoundError
        Raised **path** does not exist.

    NotADirectoryError
        Raised **path** is not a directory.

    """
    # Set the export path
    _path = path or os.getcwd()
    if not exists(_path):
        raise FileNotFoundError(f'{path} not found')
    elif not isdir(_path):
        raise NotADirectoryError(f'{path} is not a directory')

    condition = isfile if overwrite else lambda n: True
    export_dict = {k: v for k, v in EXPORT_MAPPING.items() if k in mol_format}
    if not export_dict:
        raise ValueError("No valid values found in the mol_format argument; accepted values are: "
                         "'xyz', 'pdb', 'mol' and/or 'mol2'")

    for mol in mol_list:
        mol_path = join(_path, mol.properties.name)
        if not condition:
            continue

        for ext, func in export_dict.items():
            func(mol, f'{mol_path}.{ext}')
