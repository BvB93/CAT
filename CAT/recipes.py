"""A number of recipes constructed using the CAT and Nano-CAT packages.

Examples
--------
.. code:: python

    >>> from CAT.recipes import bulk_workflow  # doctest: +SKIP
    >>> from CAT.recipes import replace_surface  # doctest: +SKIP
    >>> from CAT.recipes import dissociate_surface, row_accumulator  # doctest: +SKIP
    >>> from CAT.recipes import coordination_number  # doctest: +SKIP
    >>> from CAT.recipes import add_ligands, export_dyes, sa_scores
    ...

"""

# flake8: noqa: F403,F401

try:
    from nanoCAT import recipes as _recipes
    from nanoCAT.recipes import *

    __all__ = _recipes.__all__.copy()
    del _recipes

except ImportError as ex:
    import warnings as _warnings

    __all__ = []

    _warning = ImportWarning(str(ex))
    _warning.__cause__ = ex
    _warnings.warn(_warning)
    del _warning
    del _warnings

finally:
    from . import dye as _dye
    from .dye import *

    __all__ += _dye.__all__
    del _dye
