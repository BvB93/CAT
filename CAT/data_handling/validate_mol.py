"""A module designed for sanitizing and interpreting all molecule-related settings in the input file.

Index
-----
.. currentmodule:: CAT.data_handling.validate_mol
.. autosummary::
    santize_smiles
    validate_mol
    _parse_name_type
    _parse_mol_type

API
---
.. autofunction:: santize_smiles
.. autofunction:: validate_mol
.. autofunction:: _parse_name_type
.. autofunction:: _parse_mol_type

"""  # noqa: E501

from os.path import join, isfile, isdir, basename
from typing import Sequence, Any, Union, Optional

from rdkit import Chem
from scm.plams import Settings, Molecule
import scm.plams.interfaces.molecule.rdkit as molkit

from .validation_schemas import mol_schema
from ..utils import validate_path

__all__ = ['validate_mol', 'santize_smiles']


def santize_smiles(smiles: str) -> str:
    """Sanitize a SMILES string: turn it into a valid filename."""
    name = smiles.replace('(', '[').replace(')', ']')
    try:
        cis_trans = [item for item in smiles if item in ('/', '\\')]
        if cis_trans:
            cis_trans = [item + cis_trans[i*2+1] for i, item in enumerate(cis_trans[::2])]
            cis_trans_dict = {'//': 'trans-', '/\\': 'cis-'}
            for item in cis_trans[::-1]:
                name = cis_trans_dict[item] + name
    except KeyError:
        pass
    finally:
        return name.replace('/', '').replace('\\', '')


def validate_mol(args: Sequence[Union[Any, Settings]],
                 mol_type: str,
                 path: Optional[str] = None) -> None:
    r"""Validate the ``"input_ligands"``, ``"input_cores"`` and ``"input_qd"`` blocks in the input.

    Performs an inpalce update of **args**.

    Examples
    --------
    An example using a list of .xyz files as input

    .. code:: python

        # A list of .xyz files
        >>> args1 = ['mol1.xyz', 'mol2.xyz']
        >>> mol_type = 'input_ligands'

        >>> validate_mol(args1, mol_type)

        >>> print(args1[0], '\n', args1[1])  # doctest: +SKIP
        is_core:        False
        mol:    /path/to/my/current/working/dir/mol1.xyz
        name:   mol1
        path:   /path/to/my/current/working/dir
        type:   smiles

        is_core:        False
        mol:    /path/to/my/current/working/dir/mol2.xyz
        name:   mol2
        path:   /path/to/my/current/working/dir
        type:   smiles


    Another example using a list of .pdb-containing dictionaries as input

    .. code :: python

        # A list of Settings instances with .xyz files
        >>> print(args2)  # doctest: +SKIP
        [mol3.pdb:
             guess_bonds:       True
         mol4.pdb:
             guess_bonds:       True
        ]

        >>> mol_type = 'input_ligands'
        >>> path = '/path/to/custom/working/dir'
        >>> validate_mol(args2, mol_type, path)  # doctest: +SKIP

        >>> print(args2[0], '\n', args2[1])  # doctest: +SKIP
        guess_bonds:    True
        is_core:        True
        mol:    /path/to/custom/working/dir/mol3.pdb
        name:   mol3
        path:   /path/to/custom/working/dir
        type:   smiles

        guess_bonds:    True
        is_core:        True
        mol:    /path/to/custom/working/dir/mol4.pdb
        name:   mol4
        path:   /path/to/custom/working/dir
        type:   smiles


    Parameters
    ----------
    args : |list|_ [|plams.Settings|_]
        A list of input molecules.
        Accepts strings, PLAMS molecules and RDKit molecules.
        Additional arguments can be provided by putting above-mentioned molecules in a dictionary.

    mol_type : str
        The type of molecule.
        Accepted values are ``"input_ligands"`` and ``"input_cores"``.

    path : str
        Optional: The path to the molecule-containing directory.

    Raises
    ------
    FileNotFoundError
        Raised if **path** cannot be found.

    NotADirectoryError
        Raised if **path** is not a directory.

    ValueError
        Raised if the **mol_type** parameter is neither ``"input_cores"`` nor ``"input_ligands"``.

    SchemaError
        Raised if invalid input settings are found in while validating **args**.

    """
    # Validate arguments
    is_core = _check_core(mol_type)
    is_qd = _check_qd(mol_type)
    _path = validate_path(path)

    for i, dict_ in enumerate(args):
        if not isinstance(dict_, dict):  # No optional arguments provided
            mol = dict_
            mol_dict = Settings({'path': _path, 'is_core': is_core, 'is_qd': is_qd})
        else:  # Optional arguments have been provided: parse and validate them
            mol, mol_dict = next(iter(dict_.items()))
            mol_dict.setdefault('is_core', is_core)
            mol_dict.setdefault('is_qd', is_qd)
            mol_dict = mol_schema.validate(mol_dict)
            mol_dict.setdefault('path', _path)

        if isinstance(mol, str) and not isfile(mol):
            mol = join(mol_dict.path, mol)
        mol_dict.mol = mol

        _parse_name_type(mol_dict)
        args[i] = mol_dict


def _parse_name_type(mol_dict: Settings) -> None:
    """Set the ``"name"`` and ``"type"`` keys in **mol_dict**.

    The new values of ``"name"`` and ``"type"`` depend on the value of ``mol_dict["mol"]``.

    Parameters
    ----------
    mol_dict : |plams.Settings|_
        A Settings instance containing the ``"mol"`` key.
        ``mol_dict["mol"]`` is exp

    Raises
    ------
    TypeError
        Raised ``mol_dict["mol"]`` is an instance of neither :class:`str`, :class:`Molecule` nor
        :class:`mol`.

    """
    mol = mol_dict.mol
    if isinstance(mol, str):
        if isfile(mol):  # mol is a file
            mol_dict.type = mol.rsplit('.', 1)[-1]
            mol_dict.name = basename(mol.rsplit('.', 1)[0])
        elif isdir(mol):  # mol is a directory
            mol_dict.type = 'folder'
            mol_dict.name = basename(mol)
        else:  # mol is (probably; hopefully?) a SMILES string
            i = 1 + len(mol_dict.path) if 'path' in mol_dict else 0
            mol_dict.type = 'smiles'
            mol_dict.mol = mol[i:]
            mol_dict.name = santize_smiles(mol_dict.mol)

    elif isinstance(mol, Molecule):  # mol is an instance of plams.Molecule
        mol_dict.type = 'plams_mol'
        if not mol.properties.name:
            mol_dict.name = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)), canonical=True)
        else:
            mol_dict.name = mol.properties.name

    elif isinstance(mol, Chem.rdchem.Mol):  # mol is an instance of rdkit.Chem.Mol
        mol_dict.type = 'rdmol'
        mol_dict.name = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)

    else:
        raise TypeError(f"mol_dict['mol'] expects an instance of 'str', 'Molecule' or 'Mol'; "
                        f"observed type: '{mol.__class__.__name__}'")


def _check_core(mol_type: str) -> bool:
    """Check the **mol_type** parameter of :func:`.validate_mol`."""
    if mol_type.lower() == 'input_cores':
        return True
    elif mol_type.lower() in ('input_ligands', 'input_qd'):
        return False
    else:
        raise ValueError(f"accepted values for mol_type are 'input_cores', 'input_ligands' and "
                         f"'input_qd'; observed value: {repr(mol_type)}")


def _check_qd(mol_type: str) -> bool:
    """Parse the **mol_type** parameter of :func:`.validate_mol`."""
    if mol_type.lower() == 'input_qd':
        return True
    elif mol_type.lower() in ('input_ligands', 'input_cores'):
        return False
    else:
        raise ValueError(f"accepted values for mol_type are 'input_qd' and input_ligands; "
                         f"observed value: {repr(mol_type)}")
