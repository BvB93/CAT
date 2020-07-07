"""A module designed for optimizing the combined ligand & core.

Index
-----
.. currentmodule:: CAT.attachment.qd_opt
.. autosummary::
    init_qd_opt
    start_qd_opt
    get_job_settings
    _qd_to_db
    qd_opt

API
---
.. autofunction:: init_qd_opt
.. autofunction:: start_qd_opt
.. autofunction:: get_job_settings
.. autofunction:: _qd_to_db
.. autofunction:: qd_opt

"""

from typing import Tuple, Iterable, Optional, Type, NoReturn, Any

from scm.plams import Molecule, Settings, AMSJob
from scm.plams.core.basejob import Job

from ..jobs import job_geometry_opt  # noqa: F401
from ..workflows import WorkFlow, MOL, JOB_SETTINGS_QD_OPT, OPT
from ..mol_utils import fix_carboxyl, fix_h
from ..settings_dataframe import SettingsDataFrame
from ..data_handling.mol_to_file import mol_to_file

try:
    from nanoCAT.qd_opt_ff import qd_opt_ff
except ImportError as ex:
    def qd_opt_ff(*args: Any, ex: Exception = ex, **kwargs: Any) -> NoReturn: raise ex

__all__ = ['init_qd_opt']


def init_qd_opt(qd_df: SettingsDataFrame) -> None:
    """Initialize the ligand optimization procedure."""
    workflow = WorkFlow.from_template(qd_df, name='qd_opt')

    # Pull from the database; push unoptimized structures
    df_bool = workflow.from_db(qd_df, read_mol=True)
    idx = df_bool[OPT]
    workflow(start_qd_opt, qd_df, columns=[MOL], index=df_bool[OPT])

    # Sets a nested list
    # This cannot be done with loc is it will try to expand the list into a 2D array
    qd_df[JOB_SETTINGS_QD_OPT] = workflow.pop_job_settings(qd_df[MOL])

    # Push the optimized structures to the database
    workflow.to_db(qd_df, df_bool, columns=[JOB_SETTINGS_QD_OPT, MOL], status='optimized')

    # Export ligands to .xyz, .pdb, .mol and/or .mol format
    mol_format = qd_df.settings.optional.database.mol_format
    if mol_format:
        path = workflow.path
        mol_array = qd_df.loc[idx, MOL].values
        mol_to_file(mol_array, path, mol_format=mol_format)


def start_qd_opt(mol_list: Iterable[Molecule],
                 jobs: Tuple[Optional[Type[Job]], ...], settings: Tuple[Optional[Settings], ...],
                 use_ff: bool = False, **kwargs) -> None:
    """Loop over all molecules in **mol_list** and perform geometry optimizations."""
    for mol in mol_list:
        mol.properties.job_path = []
        qd_opt(mol, jobs, settings, use_ff=use_ff)


def qd_opt(mol: Molecule, jobs: Tuple[Optional[Type[Job]], ...],
           settings: Tuple[Optional[Settings], ...], use_ff: bool = False) -> None:
    """Perform an optimization of the quantum dot.

    Performs an inplace update of **mol**.

    Parameters
    ----------
    mol : |plams.Molecule|_
        The to-be optimized molecule.

    job_recipe : |plams.Settings|_
        A Settings instance containing all jon settings.
        Expects 4 keys: ``"job1"``, ``"job2"``, ``"s1"``, ``"s2"``.

    forcefield : bool
        If ``True``, perform the job with CP2K with a user-specified forcefield.

    """
    # Prepare the job settings
    if use_ff:
        qd_opt_ff(mol, jobs, settings)
        return None

    # Expand arguments
    job1, job2 = jobs
    s1, s2 = settings

    # Extra options for AMSJob
    if job1 is AMSJob:
        s1 = Settings(s1)
        s1.input.ams.constraints.atom = mol.properties.indices
    if job2 is AMSJob:
        s2 = Settings(s2)
        s2.input.ams.constraints.atom = mol.properties.indices

    # Run the first job and fix broken angles
    mol.job_geometry_opt(job1, s1, name='QD_opt_part1')
    fix_carboxyl(mol)
    fix_h(mol)
    mol.round_coords()

    # Run the second job
    mol.job_geometry_opt(job2, s2, name='QD_opt_part2')
    mol.round_coords()
    return None
