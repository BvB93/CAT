""" An example input file. """

import os
import yaml
import CAT.base


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = '/Users/basvanbeek/Documents/CdSe/Week_5'

# The input cores from path/core/
input_cores = yaml.load("""
-   - Cd68Se55.xyz
    - guess_bonds: False
""")

# The input ligands from path/ligand/
input_ligands = yaml.load("""
- OC(C)=O
""")

# Optional arguments: these can safely be left to their default values
argument_dict = {
    'dir_name_list': ('core', 'ligand', 'QD'),
    'dummy': 'Cl',
    'use_database': True,
    'ligand_opt': True,
    'ligand_crs': False,
    'qd_opt': False,
    'qd_int': False,
    'qd_dissociate': False,
    'split': True,
}

# Runs the script: the ligand, core and quantum dot lists are returned
qd_list, core_list, ligand_list = CAT.base.prep(input_ligands, input_cores, path, argument_dict)
