"""
Tests for the saltpy monovalent ion calculators.

Author: Irfan Alibay
Date: 2020
"""

import pytest
import math

from numpy.testing import assert_almost_equal

from saltpy import calculators


@pytest.mark.parametrize('charge, pos, neg', [
    [9, 0, 9], [0, 0, 0], [-9, 9, 0], [100000, 0, 100000]
])
def test_neutralize(charge, pos, neg):
    """Basic test for the neutralization calculator"""
    cations, anions = calculators.neutralize(charge)
    assert pos == cations
    assert neg == anions


@pytest.mark.parametrize('charge, nwat, concentration, density, pos, neg', [
    [0, 7400, 0.15, 997, 20.040, 20.040], [0, 7400, 0, 997, 0, 0],
    [0, 7400, 0.15, 10, 1998, 1998], [10, 7400, 0.15, 997, 20.040, 30.040],
    [-10, 7400, 0.15, 997, 30.040, 20.040]
])
def test_add_neutralize(charge, nwat, concentration, density, pos, neg):
    """Basic regression test for the add-then-neutralize calculator"""
    cations, anions = calculators.add_neutralize(charge, nwat, concentration,
                                                 density)
    assert_almost_equal(cations, pos, 3)
    assert_almost_equal(anions, neg, 3)


@pytest.mark.parametrize('charge, nwat, concentration, density, pos, neg', [
    [0, 7400, 0.15, 997, 20.040, 20.040], [0, 7400, 0, 997, 0, 0],
    [0, 7400, 0.15, 10, 1998, 1998], [10, 7400, 0.15, 997, 15.040, 25.040],
    [-10, 7400, 0.15, 997, 25.040, 15.040]
])
def test_split(charge, nwat, concentration, density, pos, neg):
    cations, anions = calculators.split(charge, nwat, concentration, density)
    assert_almost_equal(cations, pos, 3)
    assert_almost_equal(anions, neg, 3)


@pytest.mark.parametrize('charge, nwat, concentration, density, pos, neg', [
    [0, 7400, 0.15, 997, 20.040, 20.040], [0, 7400, 0, 997, 0, 0],
    [0, 7400, 0.15, 10, 1998, 1998], [10, 7400, 0.15, 997, 15.654, 25.654],
    [-10, 7400, 0.15, 997, 25.654, 15.654]
])
def test_sltcap(charge, nwat, concentration, density, pos, neg):
    cations, anions = calculators.sltcap(charge, nwat, concentration, density)
    assert_almost_equal(cations, pos, 3)
    assert_almost_equal(anions, neg, 3)


@pytest.mark.parametrize('conc, an, split, sltcap', [
    [0, [0, 11], [-5, 6], [0, 11]],
    [0.1, [11, 22], [6, 17], [7, 18]],
    [0.15, [16, 27], [11, 22], [12, 23]],
    [0.20, [22, 33], [17, 28], [17, 28]]
])
def test_granulysin(conc, an, split, sltcap):
    charge = 11
    nwat = 6071

    an_ions = calculators.add_neutralize(charge, nwat, conc)
    split_ions = calculators.split(charge, nwat, conc)
    sltcap_ions = calculators.sltcap(charge, nwat, conc)

    # an
    assert round(an_ions[0]) == an[0]
    assert round(an_ions[1]) == an[1]

    # split
    assert math.ceil(split_ions[0]) == split[0]
    assert math.ceil(split_ions[1]) == split[1]

    # sltcap
    assert round(sltcap_ions[0]) == sltcap[0]
    assert round(sltcap_ions[1]) == sltcap[1]


@pytest.mark.parametrize('conc, an, split, sltcap', [
    [0, [0, 20], [-10, 10], [0, 20]],
    [0.15, [20, 40], [11, 31], [12, 32]],
    [0.20, [27, 47], [17, 37], [19, 39]],
    [0.25, [33, 53], [24, 44], [25, 45]]
])
def test_jev_capsid(conc, an, split, sltcap):
    charge = 20
    nwat = 7419

    an_ions = calculators.add_neutralize(charge, nwat, conc)
    split_ions = calculators.split(charge, nwat, conc)
    sltcap_ions = calculators.sltcap(charge, nwat, conc)

    # an
    assert round(an_ions[0]) == an[0]
    assert round(an_ions[1]) == an[1]

    # split
    assert math.ceil(split_ions[0]) == split[0]
    assert math.ceil(split_ions[1]) == split[1]

    # sltcap
    assert round(sltcap_ions[0]) == sltcap[0]
    assert round(sltcap_ions[1]) == sltcap[1]
