.. image:: https://github.com/nlesc-nano/CAT/workflows/Build%20with%20Conda/badge.svg
   :target: https://github.com/nlesc-nano/CAT/actions?query=workflow%3A%22Build+with+Conda%22
.. image:: https://readthedocs.org/projects/cat/badge/?version=latest
   :target: https://cat.readthedocs.io/en/latest/
.. image:: https://zenodo.org/badge/169708827.svg
   :target: https://zenodo.org/badge/latestdoi/169708827

|

.. image:: https://img.shields.io/badge/python-3.6-blue.svg
   :target: https://docs.python.org/3.6/
.. image:: https://img.shields.io/badge/python-3.7-blue.svg
   :target: https://docs.python.org/3.7/
.. image:: https://img.shields.io/badge/python-3.8-blue.svg
   :target: https://docs.python.org/3.8/


##############################
Compound Attachment Tool 0.9.5
##############################

**CAT** is a collection of tools designed for the construction of various chemical compounds.
Further information is provided in the documentation_.

Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_.

- Create a new virtual environment, for python 3.7, using the following commands:

  - ``conda create --name CAT python``

- The virtual environment can be enabled and disabled by, respectively, typing:

  - Enable: ``conda activate CAT``

  - Disable: ``conda deactivate``


.. _dependecies:

Dependencies installation
-------------------------

Using the conda environment the following packages should be installed:

- rdkit_ : ``conda install -y --name CAT --channel conda-forge rdkit``

.. _installation:

Package installation
--------------------
Finally, install **CAT** using pip:

- **CAT**: ``pip install git+https://github.com/nlesc-nano/CAT@master --upgrade``

Now you are ready to use **CAT**.

Input files
============

Running **CAT** and can be done with the following command:
``init_cat my_settings.yaml``. The user merely has to provide a yaml_ file
with the job settings, settings which can be tweaked and altered to suit ones
purposes (see example1_). Alternatively, **CAT** can be run like a regular
python script, bypassing the command-line interface
(*i.e.* ``python input.py``, see example2_).

An extensive description of the various available settings is available in
the documentation_.


.. _yaml: https://yaml.org/
.. _documentation: https://cat.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: https://docs.anaconda.com/anaconda/install/
.. _HDF5: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QMFlows: https://github.com/SCM-NV/qmflows
.. _example1: https://github.com/BvB93/CAT/blob/master/examples/input_settings.yaml
.. _example2: https://github.com/BvB93/CAT/blob/master/examples/input.py
