.. _Optional:

Optional
========

There are a number of arguments which can be used to modify the
functionality and behavior of the quantum dot builder. Herein an
overview is provided.

Note: Inclusion of this section in the input file is not required,
assuming one is content with the default settings.

Index
~~~~~

========================================= =========================================================================================================
Option                                    Description
========================================= =========================================================================================================
:attr:`optional.database.dirname`         The name of the directory where the database will be stored.
:attr:`optional.database.read`            Attempt to read results from the database before starting calculations.
:attr:`optional.database.write`           Export results to the database.
:attr:`optional.database.overwrite`       Allow previous results in the database to be overwritten.
:attr:`optional.database.mol_format`      The file format(s) for exporting moleculair structures.
:attr:`optional.database.mongodb`         Options related to the MongoDB format.

:attr:`optional.core.dirname`             The name of the directory where all cores will be stored.
:attr:`optional.core.dummy`               Atomic number of symbol of the core dummy atoms.
:attr:`optional.core.allignment`          How the to-be attached ligands should be alligned with the core.
:attr:`optional.core.subset`              Settings related to the partial replacement of core dummy atoms.

:attr:`optional.ligand.dirname`           The name of the directory where all ligands will be stored.
:attr:`optional.ligand.optimize`          Optimize the geometry of the to-be attached ligands.
:attr:`optional.ligand.functional_groups` Manually specify SMILES strings representing functional groups.
:attr:`optional.ligand.split`             If the ligand should be attached in its entirety to the core or not.
:attr:`optional.ligand.cosmo-rs`          Perform a property calculation with COSMO-RS on the ligand.
:attr:`optional.ligand.cdft`              Perform a conceptual DFT calculation with ADF on the ligand.

:attr:`optional.qd.dirname`               The name of the directory where all quantum dots will be stored.
:attr:`optional.qd.construct_qd`          Whether or not the quantum dot should actually be constructed or not.
:attr:`optional.qd.optimize`              Optimize the quantum dot (i.e. core + all ligands).
:attr:`optional.qd.multi_ligand`          A workflow for attaching multiple non-unique ligands to a single quantum dot.
:attr:`optional.qd.bulkiness`             Calculate the :math:`V_{bulk}`, a ligand- and core-sepcific descriptor of a ligands' bulkiness.
:attr:`optional.qd.activation_strain`     Perform an activation strain analyses.
:attr:`optional.qd.dissociate`            Calculate the ligand dissociation energy.
========================================= =========================================================================================================

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    optional:
        database:
            dirname: database
            read: True
            write: True
            overwrite: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            dummy: Cl
            allignment: sphere
            subset: null

        ligand:
            dirname: ligand
            optimize: True
            functional_groups: null
            split: True
            cosmo-rs: False
            cdft: False

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            activation_strain: False
            dissociate: False
            bulkiness: False

Arguments
~~~~~~~~~

Database
--------

.. attribute:: optional.database

    All database-related settings.

    .. note::
        For :attr:`optional.database` settings to take effect the `Data-CAT <https://github.com/nlesc-nano/data-CAT>`_ package has to be installed.

    Example:

    .. code:: yaml

        optional:
            database:
                dirname: database
                read: True
                write: True
                overwrite: False
                mol_format: (pdb, xyz)
                mongodb: False

|

    .. attribute:: optional.database.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default Value** - ``"database"``

        The name of the directory where the database will be stored.

        The database directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.database.read

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Attempt to read results from the database before starting calculations.

        Before optimizing a structure, check if a geometry is available from
        previous calculations. If a match is found, use that structure and
        avoid any geometry (re-)optimizations. If one wants more control then the
        boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures will be read only for a
        specific subset.


        .. admonition:: Example

            Example #1:

            .. code:: yaml

                optional:
                    database:
                        read: (core, ligand, qd)  # This is equivalent to read: True

            Example #2:

            .. code:: yaml

                optional:
                    database:
                        read: ligand


    .. attribute:: optional.database.write

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Export results to the database.

        Previous results will **not** be overwritten unless
        :attr:`optional.database.overwrite` = ``True``. If one wants more control then
        the boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures written for a specific
        subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.overwrite

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``False``

        Allow previous results in the database to be overwritten.

        Only applicable if :attr:`optional.database.write` = ``True``.
        If one wants more control then the boolean can be substituted for
        a list of strings (*i.e.* ``"core"``, ``"ligand"`` and/or ``"qd"``), meaning
        that structures written for a specific subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.mol_format

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("pdb", "xyz")``

        The file format(s) for exporting moleculair structures.

        By default all structures are stored in the .hdf5 format as
        (partially) de-serialized .pdb files. Additional formats can be
        requested with this keyword.
        Accepted values: ``"pdb"``, ``"xyz"``, ``"mol"`` and/or ``"mol2"``.


    .. attribute:: optional.database.mongodb

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default Value** – ``False``

        Options related to the MongoDB format.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Database`:.

|

Core
----

.. attribute:: optional.core

    All settings related to the core.

    Example:

    .. code:: yaml

        optional:
            core:
                dirname: core
                dummy: Cl
                allignment: sphere
                subset: null

|

    .. attribute:: optional.core.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"core"``

        The name of the directory where all cores will be stored.

        The core directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.core.dummy

        :Parameter:     * **Type** - :class:`str` or :class:`int`
                        * **Default value** – ``17``

        Atomic number of symbol of the core dummy atoms.

        The atomic number or atomic symbol of the atoms in the core which are to be
        replaced with ligands. Alternatively, dummy atoms can be manually specified
        with the core_indices variable.


    .. attribute:: optional.core.allignment

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"sphere"``

        How the to-be attached ligands should be alligned with the core.

        Has two allowed values:

        * ``"sphere"``: Define the core vectors as those drawn from the core anchor
          atoms to the cores center.
        * ``"surface"``: Define the core vectors as those orthogonal to the cores
          surface.
          The surface is herein defined by a convex hull constructed from the core.

        Note that for a spherical core both approaches are equivalent.

        .. note::
            An example of a ``"sphere"`` (left) and ``"surface"`` (right) allignment.

            .. image:: _images/allignment.png
                :scale: 15 %
                :align: center



    .. attribute:: optional.core.subset

        :Parameter:     * **Type** - :class:`dict`, optional
                        * **Default value** – ``None``

        Settings related to the partial replacement of core dummy atoms with ligands.

        If not ``None``, has access to six further keywords,
        the first two being the most important:

        * :attr:`subset.f`
        * :attr:`subset.mode`
        * :attr:`subset.follow_edge`
        * :attr:`subset.weight`
        * :attr:`subset.randomness`
        * :attr:`subset.cluster_size`


    .. attribute:: optional.core.subset.f

        :Parameter:     * **Type** - :class:`float`

        The fraction of core dummy atoms that will actually be exchanged for ligands.

        The provided value should satisfy the following condition: :math:`0 < f \le 1`.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.core.subset.mode

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"uniform"``

        Defines how the dummy atom subset, whose size is defined by the fraction :math:`f`, will be generated.

        Accepts one of the following values:

        * ``"uniform"``: A uniform distribution; the nearest-neighbor distances between each
          successive dummy atom and all previous dummy atoms is maximized.
          can be combined with :attr:`subset.cluster_size<optional.core.subset.cluster_size>`
          to create a uniform distribution of clusters of a user-specified size.
        * ``"cluster"``: A clustered distribution; the nearest-neighbor distances between each
          successive dummy atom and all previous dummy atoms is minimized.
        * ``"random"``: A random distribution.

        It should be noted that all three methods converge towards the same set
        as :math:`f` approaches :math:`1.0`.

        If :math:`\boldsymbol{D} \in \mathbb{R}_{+}^{n,n}` is the (symmetric) distance matrix constructed
        from the dummy atom superset and :math:`\boldsymbol{a} \in \mathbb{N}^{m}` is the vector
        of indices which yields the dummy atom subset. The definition of element :math:`a_{i}`
        is defined below for the ``"uniform"`` distribution.
        All elements of :math:`\boldsymbol{a}` are furthermore constrained to be unique.

        .. math::
            :label: 1

            \DeclareMathOperator*{\argmin}{\arg\!\min}
            a_{i} = \begin{cases}
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D_{k, \hat{\imath}} \right) &
                \text{if} & i=0 \\
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{i-1} f \left( D[k, a_{\hat{\imath}}]\ \right) &
                \text{if} & i > 0
            \end{cases} \begin{matrix} & \text{with} & f(x) = e^{-x} \end{matrix}

        For the ``"cluster"`` distribution all :math:`\text{argmin}` operations
        are exchanged for :math:`\text{argmax}`.

        The old default, the p-norm with :math:`p=-2`, is equivalent to:

        .. math::
            :label: 2

            \DeclareMathOperator*{\argmax}{\arg\!\max}
            \begin{matrix}
            \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D_{k, \hat{\imath}} \right) =
            \argmax\limits_{k \in \mathbb{N}} \left( \sum_{\hat{\imath}=0}^{n} | D_{k, \hat{\imath}} |^p \right)^{1/p}
            & \text{if} & f(x) = x^{-2} \end{matrix}

        Note that as the elements of :math:`\boldsymbol{D}` were defined as positive or zero-valued real numbers;
        operating on :math:`\boldsymbol{D}` is thus equivalent to operating on its absolute.

        .. note::
            An example of a ``"uniform"``, ``"cluster"`` and ``"random"`` distribution with :math:`f=1/3`.

            .. image:: _images/distribution.png
                :scale: 15 %
                :align: center

            |
            An example of four different ``"uniform"`` distributions at :math:`f=1/16`,
            :math:`f=1/8`, :math:`f=1/4` and :math:`f=1/2`.

            .. image:: _images/distribution_p_var.png
                :scale: 20 %
                :align: center


    .. attribute:: optional.core.subset.follow_edge

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Construct the dummy atom distance matrix by following the shortest path along the
        edges of a (triangular-faced) polyhedral approximation of the core rather than the
        shortest path through space.

        Enabling this option will result in more accurate ``"uniform"`` and ``"cluster"``
        distributions at the cost of increased computational time.

        Given the matrix of Cartesian coordinates :math:`\boldsymbol{X} \in \mathbb{R}^{n, 3}`,
        the matching edge-distance matrix :math:`\boldsymbol{D}^{\text{edge}} \in \mathbb{R}_{+}^{n, n}`
        and the vector :math:`\boldsymbol{p} \in \mathbb{N}^{m}`, representing a (to-be optimized)
        path as the indices of edge-connected vertices, then element :math:`D_{i,j}^{\text{edge}}`
        is defined as following:

        .. math::
            :label: 3

            D_{i, j}^{\text{edge}} = \min_{\boldsymbol{p} \in \mathbb{N}^{m}; m \in \mathbb{N}}
            \sum_{k=0}^{m-1} || X_{p_{k},:} - X_{p_{k+1},:} ||
            \quad \text{with} \quad p_{0} = i \quad \text{and} \quad p_{m} = j

        The polyhedron edges are constructed, after projecting all vertices on the surface of a sphere,
        using Qhull's :class:`ConvexHull<scipy.spatial.ConvexHull>` algorithm
        (`The Quickhull Algorithm for Convex Hulls <https://doi.org/10.1145/235815.235821>`_).
        The quality of the constructed edges is proportional to the convexness of the core,
        more specifically: how well the vertices can be projected on a spherical surface without
        severely distorting the initial structure.
        For example, spherical, cylindrical or cuboid cores will yield reasonably edges,
        while the edges resulting from torus will be extremely poor.

        .. note::
            An example of a cores' polyhedron-representation; displaying the shortest path
            between points :math:`i` and :math:`j`.

            .. image:: _images/polyhedron.png
                :scale: 15 %
                :align: center


    .. attribute:: optional.core.subset.cluster_size

        :Parameter:     * **Type** - :class:`int` or :class:`Iterable<collections.abc.Iterable>` [:class:`int`]
                        * **Default value** – ``1``

        Allow for the creation of uniformly distributed clusters of size :math:`r`;
        should be used in conjunction with :attr:`subset.mode = "uniform"<optional.core.subset.mode>`.

        The value of :math:`r` can be either
        a single cluster size (*e.g.* :code:`cluster_size = 5`) or an iterable of various
        sizes (*e.g.* :code:`cluster_size = [2, 3, 4]`).
        In the latter case the iterable will be repeated as long as necessary.

        Compared to Eq :eq:`2` the vector of indices :math:`\boldsymbol{a} \in \mathbb{N}^{m}` is,
        for the purpose of book keeping, reshaped into the matrix
        :math:`\boldsymbol{A} \in \mathbb{N}^{q, r} \; \text{with} \; q*r = m`.
        All elements of :math:`\boldsymbol{A}` are, again, constrained to be unique.

        .. math::
            :label: 4

            \DeclareMathOperator*{\argmin}{\arg\!\min}
            A_{i,j} = \begin{cases}
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D[k, \, \hat{\imath}] \right) &
                \text{if} & i=0 & \text{and} & j=0 \\
            \argmin\limits_{k \in \mathbb{N}}
                \sum_{\hat{\imath}=0}^{i-1} \sum_{\hat{\jmath}=0}^{r} f \left( D[k, A_{\hat{\imath}, \, \hat{\jmath}}] \right) &
            \text{if} & i > 0 & \text{and} & j = 0 \\
            \argmin\limits_{k \in \mathbb{N}}
            \dfrac
                { \sum_{\hat{\imath}=0}^{i-1} \sum_{\hat{\jmath}=0}^{r} f \left( D[k, A_{\hat{\imath}, \, \hat{\jmath}}] \right) }
                { \sum_{\hat{\jmath}=0}^{j-1} f \left( D[k, A_{i, \, \hat{\jmath}}] \right) }
            &&& \text{if} & j > 0
            \end{cases}

        |

        .. note::
            An example of various cluster sizes (1, 2, 3 and 4) with :math:`f=1/4`.

            .. image:: _images/cluster_size.png
                :scale: 15 %
                :align: center

            |
            An example of clusters of varying size (:code:`cluster_size = [1, 2, 9, 1]`)
            with :math:`f=1/4`.

            .. image:: _images/cluster_size_variable.png
                :scale: 5 %
                :align: center


    .. attribute:: optional.core.subset.weight

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"numpy.exp(-x)"``

        The function :math:`f(x)` for weighting the distance.; its default value corresponds to: :math:`f(x) = e^{-x}`.

        For the old default, the p-norm with :math:`p=-2`, one can use ``weight = "x**-2"``: :math:`f(x) = x^-2`.

        Custom functions can be specified as long as they satisfy the following constraints:

        * The function must act an variable by the name of ``x``,
          a 2D array of positive and/or zero-valued floats (:math:`x \in \mathbb{R}_{+}^{n, n}`).
        * The function must take a single array as argument and return a new one.
        * The function must be able to handle values of ``numpy.nan`` and ``numpy.inf`` without
          raising exceptions.
        * The shape and data type of the output array should not change with respect to the input.

        Modules specified in the weight function will be imported when required,
        illustrated here with SciPy's :func:`expit<scipy.special.expit>`
        function: ``weight = "scipy.special.expit(x)"`` aka ``weight = "1 / (1 + numpy.exp(-x))"``

        Multi-line statements are allowed: ``weight = "a = x**2; b = 5 * a; numpy.exp(b)"``.
        The last part of the statement is assumed to be the to-be returned value
        (*i.e.* ``return numpy.exp(b)``).


    .. attribute:: optional.core.subset.randomness

        :Parameter:     * **Type** - :class:`float`, optional
                        * **Default value** – ``None``

        The probability that each new core dummy atom will be picked at random.

        Can be used in combination with ``"uniform"`` and ``"cluster"`` to introduce
        a certain degree of randomness (*i.e.* entropy).

        If not ``None``, the provided value should satisfy the following condition:
        :math:`0 \le randomness \le 1`. A value of :math:`0` is equivalent to a
        ``"uniform"`` / ``"cluster"`` distribution while :math:`1` is equivalent
        to ``"random"``.

        .. note::
            A demonstration of the ``randomness`` parameter for a ``"uniform"`` and
            ``"cluster"`` distribution at :math:`f=1/4`.

            The ``randomness`` values are (from left to right) set to :math:`0`,
            :math:`1/4`, :math:`1/2` and :math:`1`.

            .. image:: _images/randomness.png
                :scale: 13 %
                :align: center

|

Ligand
------

.. attribute:: optional.ligand

    All settings related to the ligands.

    Example:

    .. code:: yaml

        optional:
            ligand:
                dirname: ligand
                optimize: True
                functional_groups: null
                split: True
                cosmo-rs: False
                cdft: False

|

    .. attribute:: optional.ligand.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"ligand"``

        The name of the directory where all ligands will be stored.

        The ligand directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.ligand.optimize

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``True``

        Optimize the geometry of the to-be attached ligands.

        The ligand is split into one or multiple (more or less) linear fragments,
        which are subsequently optimized (RDKit UFF [1_, 2_, 3_]) and reassembled
        while checking for the optimal dihedral angle. The ligand fragments are
        biased towards more linear conformations to minimize inter-ligand
        repulsion once the ligands are attached to the core.

        After the conformation search a final (unconstrained) geometry optimization
        is performed, RDKit UFF again being the default level of theory.
        Custom job types and settings can, respectivelly, be specified with the
        ``job2`` and ``s2`` keys.

        .. note::

            .. code:: yaml

                optional:
                    ligand:
                        optimize:
                            job2: ADFJob


    .. attribute:: optional.ligand.functional_groups

        :Parameter:     * **Type** - :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** – ``None``

        Manually specify SMILES strings representing functional groups.

        For example, with :attr:`optional.ligand.functional_groups` = ``("O[H]", "[N+].[Cl-]")`` all
        ligands will be searched for the presence of hydroxides and ammonium chlorides.

        The first atom in each SMILES string (*i.e.* the "anchor") will be used for attaching the ligand
        to the core, while the last atom (assuming :attr:`optional.ligand.split` = ``True``) will be
        dissociated from the ligand and discarded.

        If not specified, the default functional groups of **CAT** are used.

        .. note::
            This argument has no value be default and will thus default to SMILES strings of the default
            functional groups supported by **CAT**.

        .. note::
            The yaml format uses ``null`` rather than ``None`` as in Python.


    .. attribute:: optional.ligand.split

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        If ``False``: The ligand is to be attached to the core in its entirety .

        =================== ==================
        Before              After
        =================== ==================
        :math:`{NR_4}^+`    :math:`{NR_4}^+`
        :math:`O_2 CR`      :math:`O_2 CR`
        :math:`HO_2 CR`     :math:`HO_2 CR`
        :math:`H_3 CO_2 CR` :math:`H_3 CO_2 CR`
        =================== ==================

        ``True``: A proton, counterion or functional group is to be removed from
        the ligand before attachment to the core.

        ========================= ==================
        Before                    After
        ========================= ==================
        :math:`Cl^- + {NR_4}^+`   :math:`{NR_4}^+`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`Na^+ + {O_2 CR}^-` :math:`{O_2 CR}^-`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`H_3 CO_2 CR`       :math:`{O_2 CR}^-`
        ========================= ==================


    .. attribute:: optional.ligand.cosmo-rs

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``


        Perform a property calculation with COSMO-RS [4_, 5_, 6_, 7_] on the ligand.

        The COSMO surfaces are by default constructed using ADF MOPAC [8_, 9_, 10_].

        The solvation energy of the ligand and its activity coefficient are
        calculated in the following solvents: acetone, acetonitrile,
        dimethyl formamide (DMF), dimethyl sulfoxide (DMSO), ethyl acetate,
        ethanol, *n*-hexane, toluene and water.


    .. attribute:: optional.ligand.cdft

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``


        Perform a conceptual DFT (CDFT) calculation with `ADF <https://www.scm.com/doc/ADF/Input/Advanced_analysis.html#conceptual-dft>`_ on the ligand.

        All global descriptors are, if installed, stored in the database.
        This includes the following properties:

        * Electronic chemical potential (mu)
        * Electronic chemical potential (mu+)
        * Electronic chemical potential (mu-)
        * Electronegativity (chi=-mu)
        * Hardness (eta)
        * Softness (S)
        * Hyperhardness (gamma)
        * Electrophilicity index (w=omega)
        * Dissocation energy (nucleofuge)
        * Dissociation energy (electrofuge)
        * Electrodonating power (w-)
        * Electroaccepting power(w+)
        * Net Electrophilicity
        * Global Dual Descriptor Deltaf+
        * Global Dual Descriptor Deltaf-

        This block can be furthermore customized with one or more of the following keys:

        * ``"keep_files"``: Whether or not to delete the ADF output afterwards.
        * ``"job1"``: The type of PLAMS Job used for running the calculation.
          The only value that should be supplied here (if any) is ``"ADFJob"``.
        * ``"s1"``: The job Settings used for running the CDFT calculation.
          Can be left blank to use the default template (:data:`nanoCAT.cdft.cdft`).

        .. admonition:: Examples

            .. code:: yaml

                optional:
                    ligand:
                        cdft: True

            .. code:: yaml

                optional:
                    ligand:
                        cdft:
                            job1: ADFJob
                            s1: ...  # Insert custom settings here

|

QD
--

.. attribute:: optional.qd

    All settings related to the quantum dots.

    Example:

    .. code:: yaml

        optional:
            qd:
                dirname: qd
                construct_qd: True
                optimize: False
                bulkiness: False
                activation_strain: False
                dissociate: False

|

    .. attribute:: optional.qd.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"qd"``

        The name of the directory where all quantum dots will be stored.

        The quantum dot directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.

    .. attribute:: optional.qd.construct_qd

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Whether or not the quantum dot should actually be constructed or not.

        Setting this to ``False`` will still construct ligands and carry out ligand workflows,
        but it will not construct the actual quantum dot itself.


    .. attribute:: optional.qd.optimize

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Optimize the quantum dot (i.e. core + all ligands) .

        By default the calculation is performed with ADF UFF [3_, 11_].
        The geometry of the core and ligand atoms directly attached to the core
        are frozen during this optimization.


    .. attribute:: optional.qd.multi_ligand

        :Parameter:     * **Type** - ``None`` or :class:`dict`
                        * **Default value** – ``None``

        A workflow for attaching multiple non-unique ligands to a single quantum dot.

        Note that this is considered a seperate workflow besides the normal ligand attachment.
        Consequently, these structures will *not* be passed to further workflows.

        See :ref:`Multi-ligand` for more details regarding the available options.

        .. note::
            An example with ``[O-]CCCC`` as main ligand and
            ``[O-]CCCCCCCCCCCCC`` & ``[O-]C`` as additional ligands.

            .. image:: _images/multi_ligand.png
                :scale: 13 %
                :align: center


    .. attribute:: optional.qd.bulkiness

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Calculate the :math:`V_{bulk}`, a ligand- and core-specific descriptor of a ligands' bulkiness.

        .. math::
            :label: 5

            V(r_{i}, h_{i}; d, h_{lim}) =
            \sum_{i=1}^{n} e^{r_{i}} (\frac{2 r_{i}}{d} - 1)^{+} (1 - \frac{h_{i}}{h_{lim}})^{+}


    .. attribute:: optional.qd.activation_strain

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Perform an activation strain analysis [12_, 13_, 14_].

        The activation strain analysis (kcal mol\ :sup:`-1`\) is performed
        on the ligands attached to the quantum dot surface with RDKit UFF [1_, 2_, 3_].

        The core is removed during this process; the analysis is thus exclusively
        focused on ligand deformation and inter-ligand interaction.
        Yields three terms:

        1.  d\ *E*\ :sub:`strain`\  : 	The energy required to deform the ligand
        from their equilibrium geometry to the geometry they adopt on the quantum
        dot surface. This term is, by definition, destabilizing. Also known as the
        preparation energy (d\ *E*\ :sub:`prep`\).

        2.  d\ *E*\ :sub:`int`\  :	The mutual interaction between all deformed
        ligands. This term is characterized by the non-covalent interaction between
        ligands (UFF Lennard-Jones potential) and, depending on the inter-ligand
        distances, can be either stabilizing or destabilizing.

        3.  d\ *E* :	The sum of d\ *E*\ :sub:`strain`\  and d\ *E*\ :sub:`int`\ .
        Accounts for both the destabilizing ligand deformation and (de-)stabilizing
        interaction between all ligands in the absence of the core.

        See :ref:`md_asa` for more details.


    .. attribute:: optional.qd.dissociate

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Calculate the ligand dissociation energy.

        Calculate the ligand dissociation energy (BDE) of ligands attached to the
        surface of the core. See :ref:`Bond Dissociation Energy` for more details.
        The calculation consists of five distinct steps:

            1.  Dissociate all combinations of |n| ligands (|Y|) and an atom from the core (|X|)
            within a radius *r* from aforementioned core atom.
            The dissociated compound has the general structure of |XYn|.

            2.  Optimize the geometry of |XYn| at the first level of theory
            (:math:`1`). Default: ADF MOPAC [1_, 2_, 3_].

            3.  Calculate the "electronic" contribution to the BDE (|dE|)
            at the first level of theory (:math:`1`): ADF MOPAC [1_, 2_, 3_].
            This step consists of single point calculations of the complete
            quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

            4.  Calculate the thermochemical contribution to the BDE (|ddG|) at the
            second level of theory (:math:`2`). Default: ADF UFF [4_, 5_]. This step
            consists of geometry optimizations and frequency analyses of the same
            compounds used for step 3.

            5.  :math:`\Delta G_{tot} = \Delta E_{1} + \Delta \Delta G_{2} = \Delta E_{1} + (\Delta G_{2} - \Delta E_{2})`.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Bond Dissociation Energy`:.



.. _1: http://www.rdkit.org
.. _2: https://github.com/rdkit/rdkit
.. _3: https://doi.org/10.1021/ja00051a040
.. _4: https://www.scm.com/doc/COSMO-RS/index.html
.. _5: https://doi.org/10.1021/j100007a062
.. _6: https://doi.org/10.1021/jp980017s
.. _7: https://doi.org/10.1139/V09-008
.. _8: https://www.scm.com/doc/MOPAC/Introduction.html
.. _9: http://openmopac.net
.. _10: https://doi.org/10.1007/s00894-012-1667-x
.. _11: https://www.scm.com/doc/UFF/index.html
.. _12: https://doi.org/10.1002/9780470125922.ch1
.. _13: https://doi.org/10.1002/wcms.1221
.. _14: https://doi.org/10.1021/acs.jpcc.5b02987

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
