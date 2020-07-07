.. _Database:

The Database Class
==================
.. currentmodule:: dataCAT

A Class designed for the storing, retrieval and updating of results.

.. image:: _images/Database.png
    :scale: 35 %
    :align: center

The methods of the Database class can be divided into three categories
accoring to their functionality:

-   Opening & closing the database - these methods serve as context managers
    for loading and unloading parts of the database from the harddrive.

    The context managers can be accessed by calling
    either :attr:`.Database.hdf5`, with the option
    of passing additional positional or keyword arguments.

    .. code:: python

        >>> from dataCAT import Database

        >>> database = Database()
        >>> with database.hdf5('r') as f:
        ...     print(type(f))
        <class 'h5py._hl.files.File'>

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    =========================  ================================
    :meth:`~Database.from_df`  :meth:`~Database.update_mongodb`
    =========================  ================================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    ======================= =
    :meth:`~Database.to_df`
    ======================= =


Index
-----

.. autosummary::

    Database.dirname
    Database.hdf5
    Database.mongodb

    Database.update_mongodb
    Database.from_df
    Database.to_df


API
---
.. autoclass:: Database
    :members:
