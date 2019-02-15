""" A module which manages all interactions with the database. """

__all__ = ['read_database', 'compare_database', 'write_database']

import os

import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit


def read_database(path, database_name='Ligand_database'):
    """
    Open the database.

    path <str>: The path to the database.
    database_name <str>: The name (including extension) of the database.

    return <pd.DataFrame>: A database of previous calculations.
    """
    path = os.path.join(path, database_name)
    if os.path.exists(path + '.json'):
        database = pd.read_json(path + '.json')
    else:
        database = pd.DataFrame()

    return database


def compare_database(plams_mol, database):
    """
    Search the database for any ligand matches.

    plams_mol <plams.Molecule>: A plams molecule.
    database <pd.DataFrame>: A database of previous calculations.

    return <plams.Molecule>, <bool>, <bool>: The (imported) ligand, if a match was found between
        input ligand and the database, and if the .pdb file of this match actually exists.
    """
    # Check if plams_mol is in database based on matches in plams.properties.name
    # Imports a molecule if a match is found
    if not database.empty and plams_mol.properties.name in list(database['Ligand_name']):
        index = list(database['Ligand_name']).index(plams_mol.properties.name)
        mol_path = list(database['Ligand_opt_pdb'])[index]
        match = True
        if os.path.exists(mol_path):
            plams_mol_new = molkit.readpdb(mol_path)
            plams_mol_new.properties = plams_mol.properties
            plams_mol = plams_mol_new
            pdb = True
        else:
            pdb = False
    else:
        match = False
        pdb = False

    return plams_mol, match, pdb


def write_database(mol_list, database, path, mol_type='ligand'):
    """
    Write the new database entries to the database.
    New entries are defined by molecules marked as mol.properties.entry = True.

    mol_list <list>[<plams.Molecule>]: A list of ligands or quantum dots.
    database <pd.DataFrame>: A database of previous calculations.
    path <str>: The path to the database.
    mol_type <str>: 'ligand' for ligands and 'qd' for quantum dots.

    return: A .json and .xlsx file.
    """
    database_entries = []
    for mol in mol_list:
        if mol.properties.entry:
            prop = mol.properties
            if mol_type == 'ligand':
                database_entries.append(
                    [prop.name,
                     prop.group,
                     mol.get_formula().split('Xx')[0],
                     os.path.join(prop.path, prop.name.split('@')[0]) + '.pdb',
                     os.path.join(prop.path, prop.name.split('@')[0]) + '.opt.pdb',
                     prop.smiles,
                     prop.energy.E_solv.Acetone,
                     prop.energy.E_solv.Acetonitrile,
                     prop.energy.E_solv.DMF,
                     prop.energy.E_solv.DMSO,
                     prop.energy.E_solv.EtOAc,
                     prop.energy.E_solv.Ethanol,
                     prop.energy.E_solv.Hexane,
                     prop.energy.E_solv.Toluene,
                     prop.energy.E_solv.Water,
                     prop.gamma_solv.Acetone,
                     prop.gamma_solv.Acetonitrile,
                     prop.gamma_solv.DMF,
                     prop.gamma_solv.DMSO,
                     prop.gamma_solv.EtOAc,
                     prop.gamma_solv.Ethanol,
                     prop.gamma_solv.Hexane,
                     prop.gamma_solv.Toluene,
                     prop.gamma_solv.Water])
            elif mol_type == 'qd':
                database_entries.append(
                    [prop.name,
                     mol.get_formula().split('Xx')[0],
                     os.path.join(prop.path, prop.name) + '.pdb',
                     os.path.join(prop.path, prop.name) + '.opt.pdb',
                     prop.energy.E,
                     prop.energy.Eint,
                     prop.energy.Estrain])

    if database_entries:
        database_entries = list(zip(*database_entries))
        if mol_type == 'ligand':
            keys = ('Ligand_name', 'Ligand_group', 'Ligand_formula', 'Ligand_pdb', 'Ligand_opt_pdb',
                    'Ligand_SMILES', 'Gsolv_Acetone', 'Gsolv_Acetonitrile', 'Gsolv_DMF',
                    'Gsolv_DMSO', 'Gsolv_EtOAc', 'Gsolv_Ethanol', 'Gsolv_Hexane', 'Gsolv_Toluene',
                    'Gsolv_Water',
                    'Gamma_Acetone', 'Gamma_Acetonitrile', 'Gamma_DMF', 'Gamma_DMSO',
                    'Gamma_EtOAc', 'Gamma_Ethanol', 'Gamma_Hexane', 'Gamma_Toluene', 'Gamma_Water')
            name = 'Ligand_database'
            sheet_name = 'Ligand'
        if mol_type == 'qd':
            keys = ('Quantum_dot_name', 'Quantum_dot_formula', 'Quantum_dot_pdb',
                    'Quantum_dot_opt_pdb', 'Quantum_dot_E', 'Quantum_dot_Eint',
                    'Quantum_dot_Estrain')
            name = 'QD_database'
            sheet_name = 'Quantum_dot'

        database_entries = pd.DataFrame(dict(zip(keys, database_entries)))

        if not database.empty:
            database = database.append(database_entries, ignore_index=True)
        else:
            database = database_entries

        path = os.path.join(path, name)
        database.to_excel(path + '.xlsx', sheet_name=sheet_name)
        database.to_json(path + '.json')
