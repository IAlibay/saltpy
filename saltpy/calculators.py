# Saltpy
# A python toolset for calculating and adding salt concentration to solvated atomistic systems
# Author: Irfan Alibay - 2020
"""
Functions to estimate monovalent salt addition --- :mod:`saltpy.calculate`
==========================================================================

References
----------
.. [Schmit2018] J. D. Schmitt, N. L. Kariyawasam, V. Needham, and P. E. Smith.
                SLTCAP: A Simple Method for Calculating the Number of Ions
                Needed for MD Simulation. J. Chem. Theory Comput.
                14(4):1823-1827, 2018.
                DOI: 10.1021/acs.jctc.7b01254
.. [Machado2020] M. R. Machado, and S. Pantano.
                 Split the Charge Difference in Two! A Rule of Thumb for
                 Adding Proper Amounts of Ions in MD Simulations.
                 J. Chem. Theory Comput. 16(3):1367-1372, 2020.
                 DOI: 10.1021/acs.jctc.9b00953
.. [Jorgensen1998] W. L. Jorgensen, and C. Jenson.
                   Temperature Dependence of TIP3P, SPC, and TIP4P Water
                   from NPT Monte Carlo Simulations: Seeking Temperatures
                   of Maximum Density. J. Comp. Chem., 19(10):1179-1186, 1998.
                   DOI: 10.1002/(SICI)1096-987X(19980730)19:10<1179::AID-JCC6>3.0.CO;2-J
"""

from __future__ import annotations

import math
from collections import namedtuple

Ions = namedtuple('ions', ['cation', 'anion'])


def neutralize(charge: float) -> tuple[int, int]:
    """
    Basic case for salt addition, only add enough monovalent ions to
    neutralize the system.

    Parameters
    ----------
    charge: float
        Total charge of solute.

    Returns
    -------
    salt: namedtuple(cations: int, anions: int)
        NamedTuple containing the cations (first entry) and anions
        (second entry) to be added to neutralize the system.
    """
    salt = Ions(charge if charge < 0 else 0,
                charge if charge > 0 else 0)
    return salt


def _water_molar_volume(nwat: int, density: float) -> float:
    """
    Helper function to calculate the molar volume of water.

    Parameters
    ----------
    nwat: int
        Number of waters in the system.
    density: float
        Density of water (g/L).

    Returns
    -------
    MolarVolume: float
        Molar volume in mol/L.
    """
    return (nwat * 18) / density


def add_neutralize(charge: int, nwat: int, concentration: float=0.15,
                   density: float=997) -> tuple[int, int]:
    """
    The Add-then-Neutralize approach to adding a given salt concentration
    to a solvated system.

    Parameters
    ----------
    charge: int
        Total charge of solute.
    nwat: int
        Number of waters in system.
    concentration: float
        Target salt molar concentration [0.15 mM].
    density: float
        Water density for molar volume calculation, default value is taken
        as the experimental water density at 298.15 K and 1 atm
        [Jorgensen1998]_ [997 g/L].

    Returns
    -------
    salt: namedtuple(cations: int, anions: int)
        NamedTuple containing the cations (first entry) and anions
        (second entry) to be added to neutralize the system.

    Notes
    -----
    This method doesn't seem to directly have a reference, however it is well
    detailed in [Schmit2018]_.
    """
    N0 = _water_molar_volume(nwat, density) * concentration
    salt = Ions(N0 + charge if charge < 0 else N0,
                N0 + charge if charge > 0 else N0)
    return salt


def split(charge: int, nwat: int, concentration: float=0.15,
          density: float=997) -> tuple[int, int]:
    """
    The SPLIT salt addition method as described in [Machado2020]_.

    Parameters
    ----------
    charge: int
        Total charge of solute.
    nwat: int
        Number of waters in system.
    concentration: float
        Target salt molar concentration [0.15 mM].
    density: float
        Water density for molar volume calculation, default value is taken
        as the experimental water density at 298.15 K and 1 atm
        [Jorgensen1998]_ [997 g/L].

    Returns
    -------
    salt: namedtuple(cations: int, anions: int)
        NamedTuple containing the cations (first entry) and anions
        (second entry) to be added to neutralize the system.

    Notes
    -----
    As detailed in [Machado2020]_, this method is an approximation of sltcap.
    See :func:`sltcap` for details on limitations.
    """
    N0 = _water_molar_volume(nwat, density) * concentration
    salt = Ions(N0 - charge/2, N0 + charge/2)
    return salt


def sltcap(charge: int, nwat: int, concentration: float=0.15,
            density: float=997) -> tuple[int, int]:
    """
    The SLTCAP method as detailed in [Schmit2018]_.

    Parameters
    ----------
    charge: int
        Total charge of solute.
    nwat: int
        Number of waters in system.
    concentration: float
        Target salt molar concentration [0.15 mM].
    density: float
        Water density for molar volume calculation, default value is taken
        as the experimental water density at 298.15 K and 1 atm
        [Jorgensen1998]_ [997 g/L].

    Returns
    -------
    salt: namedtuple(cations: int, anions: int)
        NamedTuple containing the cations (first entry) and anions
        (second entry) to be added to neutralize the system.


    Notes
    -----
    As detailed in [Schmit2018]_ and [Machado2020]_, there are several
    limitations to this method, the most notable of which are that this
    can only be used on monovalent ions. It is unclear as to how accurate
    this remains when dealing with structural ions (e.g. protein bound Ca2+).
    In this latter, case, our current approach is to treat the bound ion as a
    part of the solute (further validation is needed to verify this).
    """
    Vw = _water_molar_volume(nwat, density)
    N0 = Vw * concentration
    inner_asinh = math.asinh(charge / (2 * 1 * Vw * concentration))
    salt = Ions(N0 * math.exp(-inner_asinh),
                N0 * math.exp(inner_asinh))
    return salt

