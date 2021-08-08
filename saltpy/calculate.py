# Saltpy
# A python toolset for calculating and adding salt concentration to solvated atomistic systems
# Author: Irfan Alibay - 2020
"""
Functions to calculate salt addition from structures --- :mod:`saltpy.analyze`
==============================================================================
"""

from __future__ import annotations

from MDAnalysis.lib.mdamath import box_volume
from MDAnalysis.exceptions import NoDataError

from saltpy import estimators


def from_structure(universe: object, method: str = "SLTCAP",
                   watername: str = "HOH", concentration: float = 0.15,
                   density: float = 997):
    """
    Method for estimating the needed number of ions from an input solvated
    structure.

    Parameters
    ----------
    universe: MDAnalysis.Universe or MDAnalysis.AtomGroup
        Input solvated structure to estimate ions for. Must contain partial
        charges.
    method: str
        Method to use for estimating ions. Current options are;
        ``SLTCAP`` (default), ``SPLIT``, ``GENION``, ``ADD_NEUTRALIZE``,
        ``NEUTRALIZE``.
    watername: str
        Name of water residues ['HOH']
    concentration: float
        Target salt molar concentration [0.15 M]
    density: float
        Water density for molar volume calculations [997]. Note this is not
        used by GENION estimations.

    Returns
    -------
    salt: namedtuple(cations: int, anions: int)
        NamedTuple containing the cations (first entry) and anions
        (second entry) to be added to neutralize the system.
    """

    try:
        charge = round(universe.atoms.total_charge())
    except NoDataError as err:
        errmsg = "Input structure must contain charges"
        raise IOError(errmsg) from err

    nwat = len(universe.select_atoms(f'resname {watername}').residues)
    volume = box_volume(universe.atoms.dimensions)

    # Define a list of available methods and input arguments
    # base case kwargs set for the estimators
    kwargs_dict = {'charge': charge, 'nwat': nwat,
                   'concentration': concentration, 'density': density}
    funcs = {'SLTCAP': [estimators.sltcap, kwargs_dict],
             'SPLIT': [estimators.split, kwargs_dict],
             'ADD_NEUTRALIZE': [estimators.add_neutralize, kwargs_dict],
             'NEUTRALIZE': [estimators.neutralize, {'charge': charge}],
             'GENION': [estimators.genion, {'charge': charge, 'volume': volume,
                                            'concentration': concentration}]}

    salt = funcs[method.upper()][0](**funcs[method.upper()][1])

    return salt
