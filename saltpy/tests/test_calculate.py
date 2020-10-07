"""
Tests for the saltpy monovalent ion calculators.

Author: Irfan Alibay
Date: 2020
"""

import pytest
import math

import MDAnalysis as mda

from saltpy.calculate import from_structure
from .datafiles import TPR_6OLX, GRO_6OLX_NOIONS, TPR_2RBN, GRO_2RBN_NOIONS


@pytest.mark.parametrize('TPR, GRO, an, split, sltcap, genion', [
    [TPR_2RBN, GRO_2RBN_NOIONS, [46, 55], [42, 51], [42, 51], [49, 58]],
    [TPR_6OLX, GRO_6OLX_NOIONS, [32, 25], [29, 22], [28, 21], [36, 29]]
])
def test_calculate_regression(TPR, GRO, an, split, sltcap, genion):
    """
    Regression tests for calculate
    """
    u = mda.Universe(TPR, GRO)
    an_ions = from_structure(u, method="ADD_NEUTRALIZE", watername="SOL")
    split_ions = from_structure(u, method="SPLIT", watername="SOL")
    sltcap_ions = from_structure(u, method="SLTCAP", watername="SOL")
    genion_ions = from_structure(u, method="GENION", watername="SOL")

    # AN
    assert round(an_ions[0]) == an[0]
    assert round(an_ions[1]) == an[1]

    # SPLIT
    assert math.ceil(split_ions[0]) == split[0]
    assert math.ceil(split_ions[1]) == split[1]

    # SLTCAP
    assert round(sltcap_ions[0]) == sltcap[0]
    assert round(sltcap_ions[1]) == sltcap[1]

    # GENION
    assert round(genion_ions[0]) == genion[0]
    assert round(genion_ions[1]) == genion[1]
