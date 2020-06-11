.. _md_asa:

Ensemble-Averaged Activation Strain Analysis
============================================
.. math::
    :label: 1a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}

Herein we describe an Ensemble-Averaged extension of the
activation/strain  analysis (ASA; also known as the
distortion/interaction model), wherein the ASA is utilized
for the analyses of entire molecular dynamics trajectories.
The implementation utilizes CHARMM-style forcefields for the
calculation of all energy terms.

.. note::
    Throughout this document an overline will be used to distinguish between "normal"
    and ensemble-averaged quantities: *e.g.* :math:`E_{\text{strain}}` versus
    :math:`\overline{E}_{\text{strain}}`.

|

Strain/Distortion
-----------------
The ensemble averaged strain :math:`\Delta \overline{E}_{\text{strain}}`
represents the distortion of all ligands with respect to their equilibrium
geometry.
Given an MD trajectory with :math:`m` iterations and :math:`n` ligands per
quantum dot, the energy is averaged over all :math:`m` MD iterations and
summed over all :math:`n` ligands.

The magnitude of this term is determined by all covalent and non-covalent
intra-ligand interactions.
As this term quantifies the deviation of a ligand from its equilibrium geometry,
it is, by definition, always positive.

.. math::
    :label: 2a

    \Delta E_{\text{strain}} = E_{\text{lig-pert}} - E_{\text{lig-eq}}
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{strain}} = \frac{1}{m} \sum_{i=0}^{m} \sum_{j=0}^{n}
    E_{\text{lig-pert}}(i, j) - E_{\text{lig-eq}}


.. math::
    :label: 3a

    \Delta E_{\text{strain}} = \Delta V_{\text{bond}} + \Delta V_{\text{angle}} +
    \Delta V_{\text{Urey-Bradley}} + \Delta V_{\text{dihedral}} + \Delta V_{\text{improper}} +
    \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`E_{\text{lig-eq}}` is herein the total energy of a (single) ligand at
its equilibrium geometry, while :math:`E_{\text{lig-pert}}(i, j)` is the
total energy of the (perturbed) ligand :math:`j` at MD iteration :math:`i`.

|

Interaction
-----------
The ensemble averaged interaction :math:`\Delta \overline{E}_{\text{int}}`
represents the mutual interaction between all ligands in a molecule.
The interaction is, again, averaged over all MD iterations and
summed over all ligand-pairs.

The magnitude of this term is determined by all non-covalent inter-ligand
interactions and can be either positive (dominated by Pauli and/or
Coulombic repulsion) or negative (dominated by dispersion and/or Coulombic
attraction).

.. math::
    :label: 4a

    \Delta E_{\text{int}} = \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (j, k)
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{int}} = \frac{1}{m} \sum_{i=0}^{m}
    \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (i, j, k)

.. math::
    :label: 5a

    \Delta E_{\text{int}} = \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`\Delta E_{\text{lig-int}}(i, j, k)` represents the pair-wise
interactions between ligands :math:`j` and :math:`k` at MD iteration :math:`i`.
Double counting is avoided by ensuring that :math:`k > j`.

.. note::
    In order to avoid the substantial Coulombic repulsion between negatively charged ligands,
    its parameters are substituted with those from its neutral (*i.e.* protonated) counterpart.
    This correction is applied, exclusively, for the calculation of :math:`\Delta E_{\text{lig-int}}`.

|

Total Energy
------------
The total (ensemble-averaged) energy is the sum of
:math:`\Delta \overline{E}_{\text{strain}}` and
:math:`\Delta \overline{E}_{\text{int}}`.
Note that the energy is associated with a set of :math:`n` ligands,
*i.e.* the distortion and mutual interaction between all :math:`n` ligands.
Division by :math:`n` will thus yield the averaged energy per ligand
per MD iteration.

.. math::
    :label: 6a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}
    = \frac{1}{m} \sum_{i=0}^{m} \Delta E_{\text{strain}}(i) + \Delta E_{\text{int}}(i)

|

Examples
--------
An example input script using the ``Cd68Se55`` core and ``OC(=O)CC`` ligand.

The :attr:`activation_strain.md<optional.qd.activation_strain.md>` key enables the MD-ASA procedure;
:attr:`activation_strain.use_ff<optional.qd.activation_strain.use_ff>` ensures
that the user-specified forcefield is used during the construction of the MD trajectory.

.. code:: yaml

    path: ...

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(=O)CC

    optional:
        core:
            dummy: Cl

        ligand:
            optimize: True
            split: True

        qd:
            activation_strain:
                use_ff: True
                md: True
                job1: Cp2kJob

        forcefield:
            charge:
                keys: [input, force_eval, mm, forcefield, charge]
                Cd: 0.9768
                Se: -0.9768
                O2D2: -0.4704
                C2O3: 0.4524
            epsilon:
                unit: kjmol
                keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                Cd Cd: 0.3101
                Se Se: 0.4266
                Cd Se: 1.5225
                Cd O2D2: 1.8340
                Se O2D2: 1.6135
            sigma:
                unit: nm
                keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                Cd Cd: 0.1234
                Se Se: 0.4852
                Cd Se: 0.2940
                Cd O2D2: 0.2471
                Se O2D2: 0.3526

|

activation_strain
-----------------
.. attribute:: optional.qd.activation_strain
    :noindex:

    All settings related to the activation strain analyses.

    Example:

    .. code:: yaml

        optional:
            qd:
                activation_strain:
                    use_ff: True
                    md: True
                    iter_start: 500
                    dump_csv: False

                    el_scale14: 1.0
                    lj_scale14: 1.0

                    distance_upper_bound: "inf"
                    k: 20
                    shift_cutoff: True

                    job1: cp2kjob
                    s1: ...


            forcefield:
                ...

|

    .. attribute:: optional.qd.activation_strain.use_ff

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Utilize the parameters supplied in the :attr:`optional.forcefield` block.


    .. attribute:: optional.qd.activation_strain.md

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Perform an ensemble-averaged activation strain analysis.

        If ``True``, perform the analysis along an entire molecular dynamics trajectory.
        If ``False``, only use a single geometry instead.


    .. attribute:: optional.qd.activation_strain.iter_start

        :Parameter:     * **Type** - :class:`int`
                        * **Default value** – ``500``

        The MD iteration at which the ASA will be started.

        All preceding iteration are disgarded, treated as pre-equilibration steps.
        Note that this refers to the iteration is specified in the .xyz file.
        For example, if a geometry is written to the .xyz file very 10 iterations
        (as is the default), then ``iter_start=500`` is equivalent to
        MD iteration 5000.


    .. attribute:: optional.qd.activation_strain.dump_csv

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Dump a set of .csv files containing all potential energies gathered over the course of the MD simulation.

        For each quantum dot two files are created in the ``.../qd/asa/`` directory,
        one containing the potentials over the course of the MD simulation (``.qd.csv``) and
        for the optimized ligand (``.lig.csv``).


    .. attribute:: optional.qd.activation_strain.el_scale14

        :Parameter:     * **Type** - :class:`float`
                        * **Default value** – ``1.0``

        Scaling factor to apply to all 1,4-nonbonded electrostatic interactions.

        Serves the same purpose as the cp2k EI_SCALE14_ keyword.


    .. attribute:: optional.qd.activation_strain.lj_scale14

        :Parameter:     * **Type** - :class:`float`
                        * **Default value** – ``1.0``

        Scaling factor to apply to all 1,4-nonbonded Lennard-Jones interactions.

        Serves the same purpose as the cp2k VDW_SCALE14_ keyword.


    .. attribute:: optional.qd.activation_strain.distance_upper_bound

        :Parameter:     * **Type** - :class:`float` or :class:`str`
                        * **Default value** – ``"inf"``

        Consider only atom-pairs within this distance for calculating inter-ligand interactions.

        Units are in Angstrom.
        Using ``"inf"`` will default to the full, untruncated, distance matrix.


    .. attribute:: optional.qd.activation_strain.k

        :Parameter:     * **Type** - :class:`int`
                        * **Default value** – ``20``

        The (maximum) number of to-be considered distances per atom.

        Only relevant when :attr:`distance_upper_bound != "inf"<optional.qd.activation_strain.distance_upper_bound>`.


    .. attribute:: optional.qd.activation_strain.shift_cutoff

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Add a constant to all electrostatic and Lennard-Jones potentials such that the potential is zero at the :attr:`distance upper bound<optional.qd.activation_strain.distance_upper_bound>`.

        Serves the same purpose as the cp2k SHIFT_CUTOFF_ keyword.
        Only relevant when :attr:`distance_upper_bound != "inf"<optional.qd.activation_strain.distance_upper_bound>`.


    .. attribute:: optional.qd.activation_strain.job1

        :Parameter:     * **Type** - :class:`type` or :class:`str`
                        * **Default value** – :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass,
        used for performing the activation strain analysis.

        Should be set to :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>` if :attr:`activation_strain.md = True<optional.qd.activation_strain.md>`.


    .. attribute:: optional.qd.activation_strain.s1

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s1:
                input:
                    motion:
                        print:
                            trajectory:
                                each:
                                    md: 10
                        md:
                            ensemble: NVT
                            temperature: 300.0
                            timestep: 1.0
                            steps: 15000
                            thermostat:
                                type: CSVR
                                csvr:
                                    timecon: 1250

                    force_eval:
                        method: FIST
                        mm:
                            forcefield:
                                ei_scale14: 1.0
                                vdw_scale14: 1.0
                                ignore_missing_critical_params: ''
                                parmtype: CHM
                                parm_file_name: null
                                do_nonbonded: ''
                                shift_cutoff: .TRUE.
                                spline:
                                    emax_spline: 10e10
                                    r0_nb: 0.2
                            poisson:
                                periodic: NONE
                                ewald:
                                    ewald_type: NONE
                        subsys:
                            cell:
                                abc: '[angstrom] 100.0 100.0 100.0'
                                periodic: NONE
                            topology:
                                conn_file_format: PSF
                                conn_file_name: null
                                coord_file_format: 'OFF'
                                center_coordinates:
                                    center_point: 0.0 0.0 0.0

                    global:
                        print_level: low
                        project: cp2k
                        run_type: MD

        The job settings used for calculating the performing the ASA.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        The default settings above are specifically for the ensemble-averaged ASA
        (:attr:`activation_strain.md = True<optional.qd.activation_strain.md>`.).

.. _EI_SCALE14: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_EI_SCALE14
.. _VDW_SCALE14: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_VDW_SCALE14
.. _SHIFT_CUTOFF: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_SHIFT_CUTOFF
