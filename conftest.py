"""A pytest ``conftest.py`` file."""

import logging
from typing import Any


def pytest_configure(config: Any) -> None:
    """Flake8 is very verbose by default. Silence it.

    See Also
    --------
    https://github.com/eisensheng/pytest-catchlog/issues/59

    """
    logging.getLogger("flake8").setLevel(logging.ERROR)
    logging.getLogger("filelock").setLevel(logging.WARNING)
    logging.getLogger("CAT").setLevel(logging.WARNING)
