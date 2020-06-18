##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


0.x.0
*****
* WiP: Added an option the import pre-built quantum dots.


0.9.5
*****
* Fixed a bug where rdkit molecules were not properly converted into numpy arrays.
* Only perform `optional.ligand.optimize.job2` if the preceding UFF optimization finishes without crashing.
* Remove a ligand if its optimization fails (_i.e._ an exception is raised).
* Fixed an issue where ligand anchoring groups could not be explicitly specified
  ([docs](https://cat.readthedocs.io/en/latest/3_input_core_ligand.html#indices)).
* Improved the conformations of ligands with tri- and penta-valent pnictogen anchors.


0.9.4
*****
* Fixed an issue where certain `AMSJobs` would have duplicate keys.
* Fixed an issue where certain ligands weren't exported when constructing
  multi-ligand quantum dots (see https://github.com/nlesc-nano/CAT/issues/115).
* Re-enabled tests for Python 3.6 (see https://github.com/nlesc-nano/CAT/issues/125).


0.9.3
*****
* Switched from travis to GitHub Actions.
* Added tests using flake8, pydocstyle and doctest.


0.9.2
*****
* Updated the Database documentation (see https://github.com/nlesc-nano/data-CAT/pull/23).
* Upped the minimum Sphinx version to 2.1.
* Removed ``sphinx-autodoc-typehints``.
* Removed ``requirements.txt`` in favor of ``.readthedocs.yml``.


0.9.1
*****
* Added a new conceptual DFT (CDFT) workflow to Nano-CAT
  (https://github.com/nlesc-nano/nano-CAT/pull/57).


0.9.0
*****
* Moved a number of functions to the `nanoutils <https://github.com/nlesc-nano/Nano-Utils>`_ package.


0.8.9
*****
* Added the option to perform ligand geometry optimizations at
  arbitrary levels of theory.


0.8.8
*****
* Fixed a bug where ligands weren't properly rotated when
  using ``optional.ligand.optimize = False``.
* Added tests for ``CAT.utils.SetAttr``.


0.8.7
*****
* Replaced ``print()`` calls with ``logger.warning()`` in all dye-related functions.
* Added a recipe for calculating the synthetic accessibility score (SAS)
  as adapted from `MolGAN <https://github.com/nicola-decao/MolGAN>`_.


0.8.6
*****
* Added documentation for the new Nano-CAT ``multi_ligand_job()`` recipe.


0.8.5
*****
* Version bump.


0.8.4
*****
* Turned the ``dye`` functionality into a recipe in ``CAT.recipes``.


0.8.3
*****
* Merged all features from the ``dye`` branch into the master.
* Fixed an issue where custom forcefield settings are not properly parsed:
  https://github.com/nlesc-nano/CAT/pull/99.
* Added a try/except clause for job hashing in case rerun prevention is disabled:
  https://github.com/nlesc-nano/CAT/pull/98.
* Added new recipes to the documentation:
  https://github.com/nlesc-nano/CAT/pull/95 & https://github.com/nlesc-nano/CAT/pull/96.
* Fixed an issue where creating an object array would unpack a Molecule into Atoms:
  https://github.com/nlesc-nano/CAT/pull/94.
* Raise an Exception when failing to identify any atoms:
  https://github.com/nlesc-nano/CAT/pull/93.


0.8.2
*****
* Added the option to decorate a qd surface with more than one type of ligand.


0.8.1
*****
* Added the ``optional.core.allignment`` keyword for determining how
  ligands should be alligned with the core.
  Accepted values are ``"sphere"`` and ``"surface"``.
* https://github.com/nlesc-nano/CAT/pull/87:
  Ensure that part of the core-surface is accounted for when rotating ligands.
* https://github.com/nlesc-nano/CAT/pull/85 & https://github.com/nlesc-nano/CAT/pull/86:
  Issue a warning when atoms are too close when constructing QDs.
* https://github.com/nlesc-nano/CAT/pull/85 & https://github.com/nlesc-nano/CAT/pull/86:
  Improved warning handling.


0.8.0
*****
* Moved the ``CAT.recipes`` module to Nano-CAT.
* Moved the ``CAT.attachment.qd_opt_ff`` module to Nano-CAT.
* Created the ``CAT.workflow.key_map module`` for storing aliases
  for ``DataFrame()`` columns.
* Cleaned the modules in ``CAT.workflows``.
* Updated tests.


0.7.15
******
* Moved ``test_distribute()`` to it's own module: ``CAT.attachment.distribution_utils``.
* Added the ``brute_uniform_idx()`` for creating uniform/clustered distributions
  in a brute-force manner, *i.e.* by finding the global minimum/maximum within
  the set of all valid atom combinations.
* Generalized the ``array_combinations()`` function, it now accepts any
  array-like object and can generate combinations along any user-specified axis.
* Added the ``get_nearest_neighbors()`` function for finding the ``k``
  nearest-neighbors within a molecule.
* Added a recipe for marking a (sub-)set of surface atoms:
  ``CAT.recipes.mark_surface()``.
* Added a recipe for dissociating specific sets of surface atoms:
  ``CAT.recipes.dissociate_surface()``.
* Update to the general structure of the ``CAT.recipes`` modules.
* Multiple minor documentation adjustments.


0.7.14
******
* Changed the default value of the CP2K ``EI_SCALE14`` keyword from 0.0 to 1.0
  (*i.e.* the CHARMM forcefield default).
* Renamed the CAT ``activation_strain.scale_elstat`` keyword to ``.el_scale14``.
* Renamed the CAT ``activation_strain.scale_lj`` keyword to ``.lj_scale14``.
* Added the CAT ``activation_strain.dump_csv`` keyword for writing the raw
  potential energies to a set of .csv files.
* Added the CAT ``activation_strain.shift_cutoff`` keyword.
  Sets the value of all non-bonded potential to zero at ``activation_strain.distance_upper_bound``.
* A number of consistency improvements to the Schemas.


0.7.13
******
* Small optimization improvements to ``edge_dist()``.
* Moved a number of functions around in the CAT.utils module.
* Added the ``optional.qd.dissociate.lig_pairs`` keyword for the BDE workflow.


0.7.12
******
* Fixed a bug ``qd_opt_ff()`` where the wrong dictionary key was validated.
* Multiple updates to the CP2K MD template.
* Employ a more duck-typing based approach during the ``schema`` validation.
* Fixed a bug in the ``jobs`` module where incorrect ``Results()`` instances
  were returned.
* Multiple documentation updates.


0.7.11
******
* Updated the ``CAT.attachment.qd_opt_ff`` module in preparation for
  https://github.com/nlesc-nano/nano-CAT/pull/26.


0.7.10
******
* The function for applying distance weights during the
  subset-generation process is now configurable.
* The default distance weighting function has been changed to
  ``weight = "np.exp(-x)"``.
  The old p-norm with ``p=-2`` is still accessible via: ``weight = "x**-2"``


0.7.9
*****
* Added the option to interpolate between ``"uniform"`` / ``"cluster"`` and
  ``"random"``.
* The order of the ``p``-norm is now configurable.
* The variable representing the anchor-atom subset size has been changed
  from ``p`` to ``f``.
  ``p`` is now reserved for the order of the ``p-norm``.
* https://github.com/nlesc-nano/CAT/pull/70: Fixed an issue with the
  ``_parse_cluster_size()`` index offset.


0.7.8
*****
* It is now possible to create ``"uniform"`` distributions of clusters,
  the size of each cluster being user-specified.


0.7.7
*****
* The ``"uniform"`` and ``"cluster"`` distributions are now weighted by
  the distance rather than using a, less robust, distance truncation.


0.7.6
*****
* Added the option, when constructing core atom subsets,
  the use a distance matrix representing the shortest paths along the
  edges of a polyhedron, rather than through space.
  Enabling this option will result in more accurate ``"uniform"`` and
  ``"cluster"`` distributions at the cost of increased computational time.
* Updated and improved the ``"uniform"`` and ``"cluster"`` distributions.
* https://github.com/nlesc-nano/CAT/pull/65: Fixed a bug where ``uniform_idx()`` yielded the rolled,
  rather than unshifted, indices.
* https://github.com/nlesc-nano/CAT/pull/64: Bug fix: the subset Schema now checks for instances of
  int ``Or`` float.
* https://github.com/nlesc-nano/CAT/pull/66: Return the identity (rotation) matrix if a ``FloatingPointError`` is
  encountered during the creation of rotation matrices.
  This can occur if a ligand consists of a single atom.
* https://github.com/nlesc-nano/CAT/pull/66: Fixed a bug in the parsing of the mode parameter of ``distribute_idx()``;
  ``"uniform"`` and ``"cluster"`` will now correctly link to ``np.argmax`` and
  ``np.argmin`` instead of the other way around.


0.7.5
*****
* Added the ability to populate only a (random-ish) subset of
  core anchors with ligands.


0.7.4
*****
* The ligand rotation check is now substantially faster:
  a distance cutoff has been implemented for the construction
  of distance matrices.


0.7.3
*****
* Added an option perform an ensemble-averaged QD activation strain
  analyses in Nano-CAT_.
* Removed a number of redundant modules.
* QD optimization now properly respect the ``optional.qd.opt.use_ff`` keyword.


0.7.2
*****
* Minor tweaks to the default forcefield-related CP2K input files.
* Fixed a couple of bugs in the ligand dissociation workflow.
* Reworked the ligand dissociation procedure in Nano-CAT_.


0.7.1
*****
* Bug fix: Added a missing value to the to-be exported ASA columns.


0.7.0
*****
* Finalize the introduction of a new CAT template system (``WorkFlow()``).
* WiP: Implement an acitvation strain workflow with custom MATCH-based
  forcefields in Nano-CAT_.


0.6.5
*****
* Updated Nano-CAT to 0.2.4: https://github.com/nlesc-nano/nano-CAT/pull/20.
* Updated Data-CAT to 0.1.5: https://github.com/nlesc-nano/data-CAT/pull/17.
* Import assertions from AssertionLib_ rather than CAT_.
* Simplified to ``AsArray()`` context manager.
* Added the ``["keep_files"]`` option for quantum dot optimizations.
* Removed ``CRSJob()`` and ``CRSResults()``; import them from PLAMS_ instead.
* WiP: Introduction of a new CAT template system (``WorkFlow()``).


0.6.4
*****
* Moved the ligand bulkiness workflow from the `ligand` to the `qd` block
  in the CAT input. See `nano-CAT`_ 0.2.3.
* Updated the formula for the ligand bulkiness calculation.
  See `nano-CAT`_ 0.2.3.


0.6.3
*****
* Fixed a bug where hypervalent atoms where assigned incorrect atomic charges.


0.6.2
*****
* Added multiple improvements (and bug fixes) to the
  ligand conformation optimizer.
* Added a context manager for the `plams.Molecule.as_array()` method.
* Added an optimizer for the ligand vector.
* Updated the ligand bulkiness workflow in `nano-CAT`_ 0.2.2.


0.6.1
*****
* Added a workflow for calculating ligand bulkiness in `nano-CAT`_ 0.2.1.


0.6.0
*****
* Implemented an interface to MATCH_ (Multipurpose Atom-Typer for CHARMM)
  in Nano-CAT.
* Added a workflow for creating CP2K input files with
  the MATCH-assigned atom types & charges.
* Updated the handling of assertions, see ``CAT.assertions.assertion_manager``.


0.5.5
*****
* Lowered Python version requirement from >=3.7 to >=3.6.


0.5.4
*****
* Minor updates to the logger.
* Cleaned up CAT.jobs.py.
* ``check_sys_var()`` is now only called if an ADF-specific Job is requirest.
* Job hashes are now stored in (and retrieved from) $JN.hash files (plain text).
* Added a permanent Database_ instance to .optional.database.db.
* Parsing of functional group SMILES_ strings is now carried out during the Schema_ validation.
* Updated Data-CAT_ to 0.1.2; changed status from pre-alpha to alpha
  (see https://github.com/nlesc-nano/data-CAT/pull/13).



0.5.3
*****
* Moved Molecule to file exporting (*i.e.* .xyz and .pdb creation) from data-CAT_ to CAT_.
* Molecules can now be exported to .mol and .mol2 formats (in addition to .pdb and .xyz format).
* Increased the clarity of many exceptions (see https://github.com/nlesc-nano/CAT/issues/45).
* Updated the documentation.
* Introduced a proper logger (see https://github.com/nlesc-nano/CAT/issues/46).
* Updated data-CAT_ to 0.1.1 (https://github.com/nlesc-nano/data-CAT/pull/12) and
  nano_CAT_ to 0.1.2 (https://github.com/nlesc-nano/nano-CAT/pull/10).


0.5.2
*****
* Added more tests.
* Added a more explicit error message to ``_smiles_to_rdmol()``.


0.5.1
*****
* Documentation update.
* Updated to the ligand dissociation module in nano-CAT_ (see https://github.com/nlesc-nano/nano-CAT/issues/1).
* Added the ``keep_files`` keyword to the cosmo-rs and ligand dissociation workflows.
  Default value: ``True``.
* See https://github.com/nlesc-nano/nano-CAT/pull/9.


0.5.0
*****
* CAT_ has been split into 3 seperate packages (see https://github.com/nlesc-nano/CAT/issues/39):

  * CAT_: A collection of tools designed for the automatic construction of composite chemical compounds.
  * nano-CAT_: A collection of tools for the analysis of nanocrystals.
  * data-CAT_: A databasing framework for the Compound Attachment Tools package (CAT_).

* Docstrings have been changed into NumPy style.
* Added typehints.
* Added the CAT.SettingsDataFrame and CAT.SettingsSeries classes.
* Added more tests.
* Cleaned up all input-parsing related modules.
* Custom function groups (*i.e.* SMILES_ strings) can now be specified in the input
  under the optional.ligand.functional_groups key (see https://github.com/nlesc-nano/CAT/issues/13).


0.4.6
*****
* Added an interface between MongoDB_ and the CAT.Database_ class (see https://github.com/nlesc-nano/CAT/issues/11).


0.4.5
*****
* All raw input scripts are now stored in the structures.hdf5 file
  (see: https://github.com/nlesc-nano/CAT/issues/36).


0.4.4
*****
* Split CAT_database.py into database.py and database_functions.py.
* Unoptimized starting structures are now exported to the database.
* Added the sphinx autosummary extension.


0.4.3
*****
* Improved interaction between the database and BDE module.
* Cleaned up BDE module.
* HDF5 indices are now always sorted when itneraction with the database.


0.4.2
*****
* Numerous bug fixes.
* A couple of code-style changes.


0.4.1
*****
* COSMO-RS calculations now allow for COSMO-surface construction
  at the DFT level.


0.4.0
*****
* Introduction of the CAT.Database class.
* Central object of CAT has been changed into a dataframe of
  molecules rather than lists molecules.
* Updated a number of tests.


0.3.3
*****
* Changed qmflows template import syntax (see: https://github.com/SCM-NV/qmflows/pull/132).
* Changed yaml loader.


0.3.2
*****
* Further (minor) updates and bug fixes to the database interaction.
* Overhaul of the bond dissociation energy (BDE) module.
* Job settings are now stored in the database.


0.3.0
*****
* Massive overhaul of the CAT database interaction.
* Moved functions related to functiona group recognizition to
  CAT.attachment.ligand_anchoring.py.
* Multiple minor bug fixes.


[Unreleased]
************
* Empty Python project directory structure.


.. _AssertionLib: https://github.com/nlesc-nano/AssertionLib
.. _CAT: https://github.com/nlesc-nano/CAT
.. _CAT.Database: https://cat.readthedocs.io/en/latest/7_database.html
.. _CP2K: https://www.cp2k.org/
.. _data-CAT: https://github.com/nlesc-nano/data-CAT/
.. _Database: https://cat.readthedocs.io/en/latest/7_database.html#class-api
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _MATCH: http://brooks.chem.lsa.umich.edu/index.php?page=match&subdir=articles/resources/software
.. _MongoDB: https://www.mongodb.com/
.. _nano-CAT: https://github.com/nlesc-nano/nano-CAT/
.. _Schema: https://github.com/keleshev/schema
.. _SMILES: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
