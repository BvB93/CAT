"""Tests for the various workflows within the CAT package."""

from pathlib import Path
from typing import Optional

import yaml
import pytest
from scm.plams import Settings
from nanoutils import delete_finally, ignore_if

from CAT.base import prep

try:
    import nanoCAT  # noqa: F401
    NANOCAT_EX: Optional[ImportError] = None
except ImportError as ex:
    NANOCAT_EX = ex

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@pytest.mark.slow
@ignore_if(NANOCAT_EX)
@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_cat() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT.yaml'

    with open(yaml_path, 'r') as f:
        arg1 = Settings(yaml.load(f, Loader=yaml.FullLoader))
        arg1.path = PATH
    prep(arg1)

    with open(yaml_path, 'r') as f:
        arg2 = Settings(yaml.load(f, Loader=yaml.FullLoader))
        arg2.path = PATH
        arg2.optional.database.read = False
    prep(arg2)

    with open(yaml_path, 'r') as f:
        arg3 = Settings(yaml.load(f, Loader=yaml.FullLoader))
        arg3.path = PATH
        arg3.input_ligands += ['CSe', 'CS']
    prep(arg3)

    with open(yaml_path, 'r') as f:
        arg4 = Settings(yaml.load(f, Loader=yaml.FullLoader))
        arg4.path = PATH
        arg4.optional.database.overwrite = True
    prep(arg4)
