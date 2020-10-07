# Saltpy
# A python toolset for calculating and adding salt concentration to solvated atomistic systems
# Author: Irfan Alibay - 2020
"""
Functions to add salt concentration --- :mod:`saltpy.ionize`
============================================================
"""

import warnings
import math

import gromacs
import MDAnalysis as mda

from saltpy.calculate import from_structure


def genion(s, structure, method="SLTCAP", conc=0.15, solventname="SOL",
           **kwargs):
    """
    Wrapper around GromacsWrapper's genion call to add monovalent ions
    based on one methods implemented in :mod:`saltpy.estimators`.

    Parameters
    ----------
    s: str
        Input TPR (XDR run) file to pass to genion.
    structure: str
        Structure file matching the solvated system passed in the input TPR.
    method: str
        Estimator to use for salt addition. Can be one of: 'SLTCAP', 'SPLIT',
        'ADD-NEUTRALIZE', 'NEUTRALIZE', and 'GENION'. ['SLTCAP']
    conc: float
        Concentration of salt to be added [0.15 M].
    solventname: str
        Name of water residues in system ['SOL'].
    **kwargs
        Other GromacsWrapper genion arguments (except for `s`).
        Note, if `input` is not provided, solventname will be passed instead.

    Notes
    -----
    * If `pq` or `nq` are given non 1 & -1 values respectively, an error will
      be returned.
    * If `nn` and `np` are always ignored (mimicking the gromacs behaviour if
      `conc` is set).
    * The value of neutral is non-consequential, as all current methods will
      automatically neutralize the system.
    """

    # just in case, make sure method is in capital letters
    method = method.upper()

    # loop through kwargs, if arg mismatch is found, throw error, else
    # enforce the argument in the argument list
    for arg in [('pq', 1), ('nq', -1)]:
        if arg[0] in kwargs:
            if kwargs[arg[0]] != arg[1]:
                errmsg = (f"Argument mismatch for {arg}, value of "
                          f"{forced_inputs[arg]} was expected, got "
                          f"{kwargs[arg]} instead")
                raise TypeError(errmsg)
        else:
            kwargs[arg[0]] = arg[1]

    # check if input in kwargs, otherwise assign solventname
    if 'input' not in kwargs:
        kwargs['input'] = solventname

    # if `nn` or `np` are set, raise a warning that these will be ignored
    if ('nn' in kwargs) or ('np' in kwargs):
        wmsg = "`nn` or `np` have been passed, these will be ignored"
        warnings.warn(wmsg)

    # create the MDA universe and get the ion counts
    u = mda.Universe(s, structure)
    ions = from_structure(u, method, watername=solventname, concentration=conc)

    # add `nn` and `np` to kwargs, overwrite if necessary
    for ion in [('np', ions[0]), ('nn', ions[1])]:
        if method == 'SPLIT':
            nions = int(math.ceil(ion[1]))
        else:
            nions = int(round(ion[1]))

        if ion[0] in kwargs:
            wmsg = f"{ion[0]} passed to genion, this will be ignored"
            warnings.warn(wmsg)
        
        kwargs[ion[0]] = nions

    # run genion
    gromacs.genion(s=s, **kwargs)
