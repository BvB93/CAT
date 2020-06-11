"""A module designed for sanitizing and validating the ``["optional"]["forcefield"]`` input block.

Index
-----
.. currentmodule:: CAT.data_handling.validate_ff
.. autosummary::
    validate_ff

API
---
.. autofunction:: validate_ff

"""  # noqa: E501

import reprlib
from typing import Union, Any, Mapping, TypeVar, Optional
from collections import abc

from scm.plams import Settings, Cp2kJob

from CAT.utils import get_template

try:
    from nanoCAT.ff.cp2k_utils import set_cp2k_param
except ImportError as ex:
    NANO_CAT: Optional[ImportError] = ex
else:
    NANO_CAT: Optional[ImportError] = None

__all__ = ['validate_ff', 'update_ff_jobs']

M = TypeVar('M', bound=Mapping)


def validate_ff(s: M) -> Union[M, bool]:
    """Validate the ``["optional"]["forcefield"]`` block in the CAT input.

    Return ``False`` if an empty dictionary is provided.

    Example
    -------
    An example of a valid input block:

    .. code::

        optional:
            forcefield:
                charge:
                    keys: [input, force_eval, mm, forcefield, charge]
                    Cd: 0.9768
                    Se: -0.9768
                epsilon:
                    unit: kjmol
                    keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                    Cd Cd: 0.3101
                    Se Se: 0.4266
                    Cd Se: 1.5225
                sigma:
                    unit: nm
                    keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                    Cd Cd: 0.1234
                    Se Se: 0.4852
                    Cd Se: 0.2940

    """
    if not s:
        return False

    for k1, v1 in s.items():
        assert isinstance(k1, str), assertion_msg(k1)
        assert isinstance(v1, dict), assertion_msg(v1)

        assert 'keys' in v1
        assert isinstance(v1['keys'], abc.Sequence), assertion_msg(v1['keys'])
        v1['keys'] = tuple(v1['keys'])

        if 'unit' in v1:
            assert isinstance(v1['unit'], str), assertion_msg(v1['unit'])

        for k2, v2 in v1.items():
            if k2 in ('unit', 'keys'):
                continue
            assert isinstance(k2, str), assertion_msg(k2)
            assert isinstance(v2, (str, int, float)), assertion_msg(v2)

    return s


def assertion_msg(obj: Any) -> str:
    """Return a formatted string."""
    return f'<type {repr(obj.__class__.__name__)}: {reprlib.repr(obj)}>'


def update_ff_jobs(s: Settings) -> None:
    """Update forcefield settings."""
    if NANO_CAT is not None:
        raise NANO_CAT

    ff = Settings()
    set_cp2k_param(ff, s.optional.forcefield)

    optimize = s.optional.qd.optimize
    if optimize and optimize.use_ff:
        if optimize.job1 and str(optimize.job1) == str(Cp2kJob):
            optimize.s1 = Settings() if optimize.s1 is None else optimize.s1
            optimize.s1 += get_template('qd.yaml')['CP2K_CHARM_opt']
            optimize.s1 += ff

        if optimize.job2 and str(optimize.job1) == str(Cp2kJob):
            optimize.s2 = Settings() if optimize.s2 is None else optimize.s2
            optimize.s2 += get_template('qd.yaml')['CP2K_CHARM_opt']
            optimize.s2 += ff

    dissociate = s.optional.qd.dissociate
    if dissociate and dissociate.use_ff:
        if dissociate.job1 and str(dissociate.job1) == str(Cp2kJob):
            dissociate.s1 = Settings() if dissociate.s1 is None else dissociate.s1
            dissociate.s1 += get_template('qd.yaml')['CP2K_CHARM_opt']
            dissociate.s1 += ff

    activation_strain = s.optional.qd.activation_strain
    if activation_strain and activation_strain.use_ff:
        if activation_strain.job1 and str(activation_strain.job1) == str(Cp2kJob):
            key = 'CP2K_CHARM_singlepoint' if not activation_strain.md else 'CP2K_CHARM_md'
            activation_strain.s1 = Settings() if activation_strain.s1 is None else activation_strain.s1  # noqa
            activation_strain.s1 += get_template('qd.yaml')[key]
            activation_strain.s1 += ff
