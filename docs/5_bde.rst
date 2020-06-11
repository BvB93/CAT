.. _Bond Dissociation Energy:

Bond Dissociation Energy
========================

Calculate the bond dissociation energy (BDE) of ligands attached to the
surface of the core. The calculation consists of five distinct steps:

    1.  Dissociate all combinations of |n| ligands (|Y|, see :attr:`optional.qd.dissociate.lig_count`) a
    nd an atom from the core (|X|, see :attr:`optional.qd.dissociate.core_atom`)
    within a radius :math:`r` from aforementioned
    core atom (see :attr:`optional.qd.dissociate.lig_core_dist` and
    :attr:`optional.qd.dissociate.core_core_dist`).
    The dissociated compound has the general structure of |XYn|.

    2.  Optimize the geometry of |XYn| at the first level of theory
    (:math:`1`). Default: ADF MOPAC [1_, 2_, 3_].

    3.  Calculate the "electronic" contribution to the BDE (|dE|)
    at the first level of theory (:math:`1`): ADF MOPAC [1_, 2_, 3_].
    This step consists of single point calculations of the complete
    quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

    4.  Calculate the thermalchemical contribution to the BDE (|ddG|) at the
    second level of theory (:math:`2`). Default: ADF UFF [4_, 5_]. This step
    consists of geometry optimizations and frequency analyses of the same
    compounds used for step 3.

    5.  :math:`\Delta G_{tot} = \Delta E_{1} + \Delta \Delta G_{2} = \Delta E_{1} + (\Delta G_{2} - \Delta E_{2})`.


Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    optional:
        qd:
            dissociate:
                core_atom: Cd
                core_index: null
                lig_count: 2
                core_core_dist: 5.0  # Ångström
                lig_core_dist: 5.0  # Ångström
                lig_core_pairs: 1
                topology: {}

                keep_files: True
                job1: AMSJob
                s1: True
                job2: AMSJob
                s2: True

|

Arguments
~~~~~~~~~

.. attribute:: optional.qd.dissociate
    :noindex:

    .. code:: yaml

        optional:
            qd:
                dissociate:
                    core_atom: Cd
                    core_index: null
                    lig_count: 2
                    lig_pairs: 1
                    core_core_dist: null  # Ångström
                    lig_core_dist: 5.0  # Ångström
                    topology:
                        7: vertice
                        8: edge
                        10: face

|

    .. attribute:: optional.qd.dissociate.core_atom

        :Parameter:     * **Type** - :class:`str` or :class:`int`

        The atomic number or atomic symbol of the core atoms (:math:`X`) which are to be
        dissociated. The core atoms are dissociated in combination with :math:`n` ligands
        (:math:`Y`, see :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`).
        Yields a compound with the general formula |XYn|.

        Atomic indices can also be manually specified with :attr:`dissociate.core_index<optional.qd.dissociate.core_index>`

        If one is interested in dissociating ligands in combination with
        a molecular species (*e.g.* :math:`X = {NR_4}^+`) the atomic number (or symbol)
        can be substituted for a SMILES string represting a poly-atomic ion
        (*e.g.* tetramethyl ammonium: C[N+](C)(C)C).

        If a SMILES string is provided it must satisfy the following 2 requirements:

            1. The SMILES string *must* contain a single charged atom; unpredictable behaviour can occur otherwise.
            2. The provided structure (including its bonds) must be present in the core.

        .. warning::
            This argument has no value be default and thus *must* be provided by the user.


    .. attribute:: optional.qd.dissociate.lig_count

        :Parameter:     * **Type** - :class:`int`

        The number of ligands, :math:`n`, which is to be dissociated in combination
        with a single core atom (:math:`X`, see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`).

        Yields a compound with the general formula |XYn|.

        .. warning::
            This argument has no value be default and thus *must* be provided by the user.


    .. attribute:: optional.qd.dissociate.core_index

        :Parameter:     * **Type** - :class:`int` or :class:`list` [:class:`int`], optional
                        * **Default value** – ``None``

        Alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_core_dist>` and
        :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`.
        Manually specify the indices of all to-be dissociated atoms in the core.
        Core atoms will be dissociated in combination with the :math:`n` closest ligands.

        .. note::
            Atom numbering follows the PLAMS [1_, 2_] convention of starting from 1 rather than 0.

        .. note::
            The yaml format uses ``null`` rather than ``None`` as in Python.


    .. attribute:: optional.qd.dissociate.core_core_dist

        :Parameter:     * **Type** - :class:`float` or :class:`int`, optional
                        * **Default value** – ``None``

        The maximum to be considered distance (Ångström) between atoms in
        :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`.
        Used for determining the topology of the core atom

        (see :attr:`dissociate.topology<optional.qd.dissociate.topology>`) and whether it is exposed to the
        surface of the core or not. It is recommended to use a radius which
        encapsulates a single (complete) shell of neighbours.

        If not specified (or equal to ``0.0``) **CAT** will attempt to guess a suitable value
        based on the cores' radial distribution function.


    .. attribute:: optional.qd.dissociate.lig_core_dist

        :Parameter:     * **Type** - :class:`float` or :class:`int`, optional
                        * **Default value** – ``None``

        Dissociate all combinations of a single core atom (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`)
        and the :math:`n` closests ligands within a user-specified radius.

        Serves as an alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_pairs>`,
        which removes a set number of combinations rather than everything withing a certain radius.

        The number of ligands dissociated in combination with a single core atom is controlled by
        :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`.

        .. image:: _images/BDE_XY2.png
            :scale: 25 %
            :align: center

|


    .. attribute:: optional.qd.dissociate.lig_pairs

        :Parameter:     * **Type** - :class:`int`, optional
                        * **Default value** – ``None``

        Dissociate a user-specified number of combinations of a single core atom (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`)
        and the :math:`n` closests ligands.

        Serves as an alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_core_dist>`,
        removing a preset number of (closest) pairs rather than all combinations within a certain radius.

        The number of ligands dissociated in combination with a single core atom is controlled by
        :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`.


    .. attribute:: optional.qd.dissociate.topology

        :Parameter:     * **Type** - :class:`dict`
                        * **Default value** – ``{}``

        A dictionary which translates the number neighbouring core atoms
        (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>` and
        :attr:`dissociate.core_core_dist<optional.qd.dissociate.core_core_dist>`)
        into a topology. Keys represent the number of neighbours, values represent
        the matching topology.

        .. admonition:: Example

            Given a :attr:`dissociate.core_core_dist<optional.qd.dissociate.core_core_dist>` of ``5.0`` Ångström,
            the following options can be interpreted as following:

            .. code:: yaml

                optional:
                    qd:
                        dissociate:
                            7: vertice
                            8: edge
                            10: face

            Core atoms with ``7`` other neighbouring core atoms (within a radius of ``5.0`` Ångström)
            are marked as ``"vertice"``, the ones with ``8`` neighbours are marked as ``"edge"``
            and the ones with ``10`` neighbours as ``"face"``.

|

Arguments - Job Customization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. attribute:: optional.qd.dissociate
    :noindex:

    .. code:: yaml

        optional:
            qd:
                dissociate:
                    keep_files: True
                    job1: AMSJob
                    s1: True
                    job2: AMSJob
                    s2: True

|

    .. attribute:: optional.qd.dissociate.keep_files

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Whether to keep or delete all BDE files after all calculations are finished.


    .. attribute:: optional.qd.dissociate.job1

        :Parameter:     * **Type** - :class:`type`, :class:`str` or :class:`bool`
                        * **Default value** – :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass, used for calculating the
        "electronic" component (|dE_lvl1|) of the bond dissociation energy.
        Involves single point calculations.

        Alternatively, an alias can be provided for a specific
        job type (see :ref:`Type Aliases`).

        Setting it to ``True`` will default to :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`,
        while ``False`` is equivalent to :attr:`optional.qd.dissociate` = ``False``.


    .. attribute:: optional.qd.dissociate.s1

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s1:
                input:
                    mopac:
                        model: PM7
                    ams:
                        system:
                            charge: 0

        The job settings used for calculating the "electronic" component
        (|dE_lvl1|) of the bond dissociation energy.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        Setting it to ``True`` will default to the ``["MOPAC"]`` block in
        CAT/data/templates/qd.yaml_, while ``False`` is equivalent to
        :attr:`optional.qd.dissociate` = ``False``.


    .. attribute:: optional.qd.dissociate.job2

        :Parameter:     * **Type** - :class:`type`, :class:`str` or :class:`bool`
                        * **Default value** – :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass, used for calculating the
        thermal component (|ddG_lvl2|) of the bond dissociation energy.
        Involves a geometry reoptimizations and frequency analyses.

        Alternatively, an alias can be provided for a specific
        job type (see :ref:`Type Aliases`).


        Setting it to ``True`` will default to :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`,
        while ``False`` will skip the thermochemical analysis completely.


    .. attribute:: optional.qd.dissociate.s2

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s2:
                input:
                    uff:
                        library: uff
                    ams:
                        system:
                            charge: 0
                            bondorders:
                                _1: null

        The job settings used for calculating the thermal component (|ddG_lvl2|)
        of the bond dissociation energy.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        Setting it to ``True`` will default to the the *MOPAC* block in
        CAT/data/templates/qd.yaml_, while ``False`` will skip the
        thermochemical analysis completely.

Index
-----
.. currentmodule:: nanoCAT.bde.dissociate_xyn
.. autosummary::
    dissociate_ligand
    MolDissociater
    MolDissociater.remove_bulk
    MolDissociater.assign_topology
    MolDissociater.get_pairs_closest
    MolDissociater.get_pairs_distance
    MolDissociater.combinations
    MolDissociater.__call__

API
---
.. autofunction:: dissociate_ligand
.. autoclass:: MolDissociater
.. automethod:: MolDissociater.remove_bulk
.. automethod:: MolDissociater.assign_topology
.. automethod:: MolDissociater.get_pairs_closest
.. automethod:: MolDissociater.get_pairs_distance
.. automethod:: MolDissociater.combinations
.. automethod:: MolDissociater.__call__


.. _1: https://www.scm.com/doc/MOPAC/Introduction.html
.. _2: http://openmopac.net
.. _3: https://doi.org/10.1007/s00894-012-1667-x
.. _4: https://doi.org/10.1021/ja00051a040
.. _5: https://www.scm.com/doc/UFF/index.html
.. _qd.yaml: https://github.com/BvB93/CAT/blob/master/CAT/data/templates/qd.yaml

.. |dE| replace:: :math:`\Delta E`
.. |dE_lvl1| replace:: :math:`\Delta E_{1}`
.. |dE_lvl2| replace:: :math:`\Delta E_{2}`
.. |dG| replace:: :math:`\Delta G_{tot}`
.. |dG_lvl2| replace:: :math:`\Delta G_{2}`
.. |ddG| replace:: :math:`\Delta \Delta G`
.. |ddG_lvl2| replace:: :math:`\Delta \Delta G_{2}`
.. |XYn| replace:: :math:`XY_{n}`
.. |Yn| replace:: :math:`Y_{n}`
.. |n| replace:: :math:`{n}`
.. |X| replace:: :math:`X`
.. |Y| replace:: :math:`Y`
