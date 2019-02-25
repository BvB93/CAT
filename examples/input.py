""" An example input file. """

import yaml
import CAT
from scm.plams.core.settings import Settings


yml = '/Users/basvanbeek/Documents/GitHub/CAT/examples/input_settings.yaml'
with open(yml, 'r') as file:
    arg = Settings(yaml.load(file))

qd_list, core_list, ligand_list = CAT.base.prep(arg)
