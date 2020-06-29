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
    either :attr:`.Database.csv_lig`, :attr:`.Database.csv_qd`,
    or :attr:`.Database.hdf5`, with the option
    of passing additional positional or keyword arguments.

    .. code:: python

        >>> from dataCAT import Database

        >>> database = Database()
        >>> with database.csv_lig(write=False) as db:
        >>>     print(repr(db))
        DFProxy(ndframe=<pandas.core.frame.DataFrame at 0x7ff8e958ce80>)

        >>> with database.hdf5('r') as db:
        >>>     print(type(db))
        <class 'h5py._hl.files.File'>

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    ============================  =============================  ================================
    :meth:`~Database.update_csv`  :meth:`~Database.update_hdf5`  :meth:`~Database.update_mongodb`
    ============================  =============================  ================================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    ==========================  ===========================
    :meth:`~Database.from_csv`  :meth:`~Database.from_hdf5`
    ==========================  ===========================


Index
-----

.. autosummary::

    Database.dirname
    Database.csv_lig
    Database.csv_qd
    Database.hdf5
    Database.mongodb

    Database.update_mongodb
    Database.update_csv
    Database.update_hdf5
    Database.from_csv
    Database.from_hdf5

    DFProxy
    OpenLig
    OpenQD


API
---
.. autoclass:: Database
    :members:

.. autoclass:: DFProxy

.. autoclass:: OpenLig
    :members:
    :inherited-members:

.. autoclass:: OpenQD
    :members:
    :inherited-members:
