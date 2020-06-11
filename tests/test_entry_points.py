"""Tests for :mod:`CAT.data_handling.entry_points`."""

from pathlib import Path
from assertionlib import assertion
from nanoutils import delete_finally

from CAT.data_handling.entry_points import main

PATH = Path('tests') / 'test_files'
LIG_PATH = PATH / 'ligand'
QD_PATH = PATH / 'qd'
DB_PATH = PATH / 'database'


@delete_finally(LIG_PATH, QD_PATH, DB_PATH)
def test_main() -> None:
    """Test :func:`CAT.data_handling.entry_points.main`."""
    filename = str(PATH / 'input2.yaml')
    main([filename])

    assertion.assert_(main, [f'{filename}bob'], exception=FileNotFoundError)
