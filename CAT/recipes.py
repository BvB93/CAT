"""
CAT.recipes
===========

A number of recipes constructed using the CAT and Nano-CAT packages.

Examples
--------
.. code:: python

    >>> from CAT.recipes import bulk_workflow
    >>> from CAT.recipes import replace_surface
    >>> from CAT.recipes import dissociate_surface, row_accumulator
    >>> from CAT.recipes import coordination_number
    >>> from CAT.recipes import add_ligands, export_dyes, sa_scores
    >>> from CAT.recipes import mutli_ligand_job
    ...

"""

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
