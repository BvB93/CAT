#!/usr/bin/env python
"""Entry points for CAT.

Index
-----
.. currentmodule:: CAT.data_handling.input_parser
.. autosummary::
    extract_args
    main

API
---
.. autofunction:: extract_args
.. autofunction:: main

"""

import argparse
from os import getcwd
from os.path import join, exists
from typing import Optional, List

import yaml

from scm.plams import Settings

from CAT import base
from CAT.logger import logger


def extract_args(args: Optional[List[str]] = None) -> Settings:
    """Extract and return all arguments."""
    input_file = args.YAML[0]
    if exists(input_file):
        pass
    elif exists(join(getcwd(), input_file)):
        input_file = join(getcwd(), input_file)
    else:
        input_file2 = join(getcwd(), input_file)
        raise FileNotFoundError(f'No file found at {input_file} or {input_file2}')

    with open(input_file, 'r') as file:
        return Settings(yaml.load(file, Loader=yaml.FullLoader))


def main(args: Optional[List[str]] = None) -> None:
    """Launch CAT via the command line."""
    parser = argparse.ArgumentParser(
        prog='CAT',
        usage='init_cat my_settings_file.yaml',
        description=('Description: This script initalizes the Compound Attachment Tool (CAT).')
    )

    parser.add_argument(
        'YAML', nargs=1, type=str, metavar='input.yaml',
        help='Required: A .yaml file with the settings for CAT'
    )

    try:
        args = parser.parse_args(args)
        base.prep(extract_args(args), return_mol=False)
    except Exception as ex:
        logger.critical(f'{ex.__class__.__name__}: {ex}')
        raise ex
