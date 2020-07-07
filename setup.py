#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit CAT/__version__.py
version = {}
with open(os.path.join(here, 'CAT', '__version__.py')) as f:
    exec(f.read(), version)
__version__: str = version['__version__']

with open('README.rst', encoding='utf-8') as readme_file:
    readme = readme_file.read()

# Check if rdkit is manually installed (as it is not available via pypi)
try:
    import rdkit
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(
        "'CAT' requires the 'rdkit' package: https://anaconda.org/conda-forge/rdkit"
    ) from None

# Check the dataCAT version
try:
    import dataCAT
except ModuleNotFoundError:
    pass
else:
    version_info = getattr(dataCAT, 'version_info', (-1, -1, -1))
    if version_info < (0, 6):
        raise ImportError(f"'CAT' {__version__} requires 'dataCAT' >= 0.6.0; "
                          f"observed version: {dataCAT.__version__}")

docs_require = [
    'sphinx>=2.4,<3.1',
    'sphinx_rtd_theme',
    'data-CAT@git+https://github.com/nlesc-nano/data-CAT@devel'
]

tests_require = [
    'pytest>=5.4.0',
    'pytest-cov',
    'flake8>=3.8.0',
    'pyflakes>=2.1.1',
    'pytest-flake8>=1.0.6',
    'pydocstyle>=5.0.0',
    'pytest-pydocstyle>=2.1',
    'pytest-mock'
]

if sys.version_info[1] > 6:
    docs_require.append('nano-CAT@git+https://github.com/nlesc-nano/nano-CAT@devel')
    docs_require.append('auto-FOX@git+https://github.com/nlesc-nano/auto-FOX@devel')
tests_require += docs_require

setup(
    name='CAT',
    version=__version__,
    description=('A collection of tools designed for the automatic '
                 'construction of chemical compounds.'),
    long_description=f'{readme}\n\n',
    long_description_content_type='text/x-rst',
    author=['Bas van Beek', 'Jelena Belic'],
    author_email='b.f.van.beek@vu.nl',
    url='https://github.com/nlesc-nano/CAT',
    packages=[
        'CAT',
        'CAT.attachment',
        'CAT.data',
        'CAT.dye',
        'CAT.data.coskf',
        'CAT.data.templates',
        'CAT.data.CORES',
        'CAT.data.LIGANDS',
        'CAT.data_handling',
        'CAT.workflows'
    ],
    package_dir={'CAT': 'CAT'},
    package_data={
        'CAT': [
            'data/templates/*yaml',
            'data/coskf/*coskf',
            'data/CORES/*xyz',
            'data/LIGANDS/*xyz',
            'workflows/workflow_yaml.yaml',
            'workflows/*.pyi',
            'py.typed'
        ]
    },
    entry_points={
        'console_scripts': ['init_cat=CAT.data_handling.entry_points:main']
    },
    include_package_data=True,
    license='GNU Lesser General Public License v3 or later',
    zip_safe=False,
    keywords=[
        'quantum-mechanics',
        'molecular-mechanics',
        'science',
        'chemistry',
        'python-3',
        'python-3-6',
        'python-3-7',
        'python-3-8',
        'automation',
        'scientific-workflows'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Chemistry'
        'Typing :: Typed'
    ],
    test_suite='tests',
    python_requires='>=3.6',
    install_requires=[
        'h5py',
        'Nano-Utils>=0.4.3',
        'numpy',
        'scipy',
        'pandas<1.0.0',
        'pyyaml>=5.1',
        'schema',
        'AssertionLib>=2.2.3',
        'plams@git+https://github.com/SCM-NV/PLAMS@a5696ce62c09153a9fa67b2b03a750913e1d0924',
        'qmflows@git+https://github.com/SCM-NV/qmflows@master',
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=tests_require,
    extras_require={
        'test': tests_require,
        'doc': docs_require
    }
)
