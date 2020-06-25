.. _PDBContainer:

The PDBContainer Class
======================
.. currentmodule:: dataCAT

An immutable class for holding array-like representions of a set of .pdb files.

Index
-----

.. autosummary::

    PDBContainer.atoms
    PDBContainer.bonds
    PDBContainer.atom_count
    PDBContainer.bond_count
    PDBContainer.index

    PDBContainer.__init__
    PDBContainer.__getitem__
    PDBContainer.__len__
    PDBContainer.keys
    PDBContainer.values
    PDBContainer.items

    PDBContainer.from_molecules
    PDBContainer.to_molecules
    PDBContainer.create_hdf5_group
    PDBContainer.validate_hdf5
    PDBContainer.from_hdf5
    PDBContainer.to_hdf5


API
---
.. autoclass:: PDBContainer
    :members: atoms, bonds, atom_count, bond_count, index, __init__

.. automethod:: PDBContainer.__getitem__
.. automethod:: PDBContainer.__len__
.. automethod:: PDBContainer.keys
.. automethod:: PDBContainer.values
.. automethod:: PDBContainer.items
.. automethod:: PDBContainer.from_molecules
.. automethod:: PDBContainer.to_molecules
.. automethod:: PDBContainer.create_hdf5_group
.. automethod:: PDBContainer.validate_hdf5
.. automethod:: PDBContainer.from_hdf5
.. automethod:: PDBContainer.to_hdf5
