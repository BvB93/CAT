"""Tests for :mod:`CAT.thermo_chem`."""

from os.path import join

import numpy as np

from scm.plams import (Molecule, Units)
from assertionlib import assertion

from CAT.thermo_chem import (get_entropy, get_thermo)

PATH = join('tests', 'test_files')
MOL = Molecule(join(PATH, 'Methanol.xyz'))  # Methanol; BP86/QZ4P
FREQ = np.load(join(PATH, 'freq.npy'))  # cm**-1


def test_get_entropy() -> None:
    """Tests for :func:`CAT.thermo_chem.get_entropy`."""
    ref1 = np.array([143.78052972, 81.90558458, 2308.88449109])
    s1 = get_entropy(MOL, FREQ)
    np.testing.assert_allclose(ref1, s1)

    # Test with a different temperature
    ref2 = np.array([149.8889032, 85.57060867, 2338.20467688])
    s2 = get_entropy(MOL, FREQ, T=400)
    np.testing.assert_allclose(ref2, s2)


def test_get_thermo() -> None:
    """Tests for :func:`CAT.thermo_chem.get_thermo`."""
    ref1 = {'E': 0.0,
            'U': 36.876336257248056,
            'H': 37.46882030779777,
            'S': 0.07643865635175046,
            'G': 14.678634916523373}
    thermo1 = get_thermo(MOL, FREQ)
    for k, v in thermo1.items():
        i, j = round(v, 8), round(ref1[k], 8)
        assertion.eq(i, j)

    # Test with E != 0.0
    ref2 = {'E': -692.08,
            'U': -655.203663742752,
            'H': -654.6111796922023,
            'S': 0.07643865635175046,
            'G': -677.4013650834768}
    thermo2 = get_thermo(MOL, FREQ, -692.08)
    for k, v in thermo2.items():
        i, j = round(v, 8), round(ref2[k], 8)
        assertion.eq(i, j)

    # Test with a different unit (au)
    thermo3 = get_thermo(MOL, FREQ, unit='au')
    for k, v in thermo3.items():
        i, j = round(v, 8), round(Units.convert(ref1[k], 'kcal/mol', 'au'), 8)
        assertion.eq(i, j)

    # Test with a different temperature
    ref4 = {'E': 0.0,
            'U': 39.09512888121028,
            'H': 39.89000937834221,
            'S': 0.08340885418719352,
            'G': 6.526467703464801}
    thermo4 = get_thermo(MOL, FREQ, T=400)
    for k, v in thermo4.items():
        i, j = round(v, 8), round(ref4[k], 8)
        assertion.eq(i, j)

    # Test when exporting a single quantity
    ref5 = 14.678634916523373
    g = get_thermo(MOL, FREQ, export='G')
    i, j = round(g, 8), round(ref5, 8)
    assertion.eq(i, j)
