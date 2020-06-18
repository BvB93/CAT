"""Tests for the various workflows within the CAT package."""

from pathlib import Path

import yaml
from scm.plams import Settings
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.base import prep
from CAT.workflows import MOL

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_cat() -> None:
    """Tests for the CAT package."""
    yaml_path = PATH / 'CAT_indices.yaml'
    with open(yaml_path, 'r') as f:
        arg = Settings(yaml.load(f, Loader=yaml.FullLoader))

    arg.path = PATH
    *_, ligand_df = prep(arg)

    ref = {('CCCOCCC(=O)O', 'C5'),
           ('CCCNCCC(=O)O', 'C5'),
           ('CCCCCCC(=O)[O-]', 'O8')}
    assertion.eq(set(ligand_df.index), ref)

    for (name, i), mol in ligand_df[MOL].items():
        symbol = i[0]
        idx = int(i[1])
        atom = mol[idx]

        assertion.is_(mol.properties.dummies, atom, message=name)
        assertion.eq(symbol, atom.symbol, message=name)
