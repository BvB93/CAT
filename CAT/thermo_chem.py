"""A module related to calculating thermochemical properties.

Index
-----
.. currentmodule:: CAT.thermo_chem
.. autosummary::
    get_entropy
    get_thermo

API
---
.. autofunction:: get_entropy
.. autofunction:: get_thermo

"""

from typing import (Sequence, Union, Dict)

import numpy as np

from scm.plams import Molecule
from scm.plams.tools.units import Units

__all__ = ['get_thermo', 'get_entropy']

# flake8: noqa: N806,N803


def get_entropy(mol: Molecule,
                freqs: np.ndarray,
                T: float = 298.15) -> np.ndarray:
    """Calculate the translational, vibrational and rotational entropy.

    All units and constants are in SI units.

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    freqs : |np.ndarray|_ [|np.float64|_]
        An iterable consisting of vibrational frequencies in units of cm**-1.

    T : float
        The temperature in Kelvin.

    Returns
    -------
    |np.ndarray|_ [|np.float64|_]:
        An array with translational, rotational and vibrational contributions to the entropy,
        ordered in that specific manner.
        Units are in J/mol.

    """
    # Define constants
    kT = 1.380648 * 10**-23 * T  # Boltzmann constant * temperature
    h = 6.6260701 * 10**-34  # Planck constant
    hv_kT = (h * np.asarray(freqs)) / kT  # (Planck * frequencies) / (Boltzmann * temperature)
    R = 8.31445  # Gas constant
    V_Na = ((R * T) / 10**5) / Units.constants['NA']  # Volume(1 mol ideal gas) / Avogadro's number
    pi = np.pi

    # Extract atomic masses and Cartesian coordinates
    m = np.array([at.mass for at in mol]) * 1.6605390 * 10**-27
    x, y, z = mol.as_array().T * 10**-10

    # Calculate the rotational partition function: q_rot
    inertia = np.array([
        [sum(m*(y**2 + z**2)), -sum(m*x*y), -sum(m*x*z)],
        [-sum(m*x*y), sum(m*(x**2 + z**2)), -sum(m*y*z)],
        [-sum(m*x*z), -sum(m*y*z), sum(m*(x**2 + y**2))]
    ])
    inertia_product = np.product(np.linalg.eig(inertia)[0])
    q_rot = pi**0.5 * ((8 * pi**2 * kT) / h**2)**1.5 * inertia_product**0.5

    # Calculate the translational, rotational and vibrational entropy (divided by R)
    S_trans = 1.5 + np.log(V_Na * ((2 * pi * sum(m) * kT) / h**2)**1.5)
    S_rot = 1.5 + np.log(q_rot)
    with np.errstate(divide='ignore', invalid='ignore'):
        S_vib_left = hv_kT / np.expm1(hv_kT)
        S_vib_left[np.isnan(S_vib_left)] = 0.0
        S_vib_right = np.log(1 - np.exp(-hv_kT))
        S_vib_right[S_vib_right == -np.inf] = 0.0
    S_vib = sum(S_vib_left - S_vib_right)

    return R * np.array([S_trans, S_rot, S_vib])


def get_thermo(mol: Molecule,
               freqs: Sequence[float],
               E: float = 0.0,
               T: float = 298.15,
               export: Sequence[str] = ('E', 'U', 'H', 'S', 'G'),
               unit: str = 'kcal/mol') -> Union[float, Dict[str, float]]:
    """Extract and return Gibbs free energies, entropies and/or enthalpies from an AMS KF file.

    All vibrational frequencies smaller than 100 cm**-1 are set to 100 cm**-1.

    .. _plams.Units: https://www.scm.com/doc/plams/components/utils.html#scm.plams.tools.units.Units

    Parameters
    ----------
    mol : |plams.Molecule|_
        A PLAMS molecule.

    freqs : |np.ndarray|_ [|np.float64|_]
        An iterable consisting of vibrational frequencies in units of cm**-1.

    E : float
        The eletronic energy in kcal/mol.
        Defaults to 0.0 kcal/mol.

    T : float
        The temperature in Kelvin.

    export : |tuple|_ [|str|_]
        An iterable containing strings of the to be exported energies:

        * ``'E'``: Electronic energy (see the **E** parameter)
        * ``'U'``: Interal energy (:math:`E + U_{nuc}`)
        * ``'H'``: Enthalpy (:math:`U + pV`)
        * ``'S'``: Entropy
        * ``'G'``: Gibbs free energy (:math:`H - T*S`)

    unit : str
        The unit of the to be returned energies.
        See plams.Units_ for more details and an overview of available energy units.

    Returns
    -------
    |float|_ or |dict|_ [|str|_, |float|_]:
        An energy or dictionary of energies.
        Keys

    """
    # Get frequencies; set all frequencies smaller than 100 cm**-1 to 100 cm**-1
    freqs = np.array(freqs)
    freqs[freqs < 100] = 100
    freqs *= 100 * Units.constants['c']

    # hv_kT = (Planck constant * frequencies) / (Boltzmann constant * temperature)
    hv_kT = (6.6260701 * 10**-34 * freqs) / (1.380648 * 10**-23 * T)
    RT = 8.31445 * T  # Gas constant * temperature

    # Extract and/or calculate the various energies
    E = E * Units.conversion_ratio('kcal/mol', 'kj/mol') * 1000
    U = E + RT * (3.0 + sum(0.5 * hv_kT + hv_kT / np.expm1(hv_kT)))
    H = U + RT
    S = sum(get_entropy(mol, freqs, T=T))
    G = H - T * S

    ret = {'E': E, 'U': U, 'H': H, 'S': S, 'G': G}

    if len(export) == 1:
        return Units.convert(ret[export[0]], 'kj/mol', unit) / 1000
    return {i: Units.convert(ret[i], 'kj/mol', unit) / 1000 for i in ret if i in export}
