"""
Tests for the saltpy monovalent ion calculators.

Author: Irfan Alibay
Date: 2020
"""

import pytest
import math

import MDAnalysis as mda

from saltpy.ionize import genion
from .datafiles import (TPR_6OLX, GRO_6OLX_NOIONS,
                        TPR_2RBN, GRO_2RBN_NOIONS,)


@pytest.mark.parametrize('method, pos, neg', [
    ['ADD_NEUTRALIZE', 46, 55],
    ['SPLIT', 42, 51],
    ['SLTCAP', 42, 51],
    ['GENION', 49, 58]
])
def test_2rbn(tmpdir, method, pos, neg):
    with tmpdir.as_cwd():
        genion(s=TPR_2RBN, structure=GRO_2RBN_NOIONS, o='test.gro',
               method=method)
        test = mda.Universe('test.gro')
        assert len(test.select_atoms('resname NA')) == pos
        assert len(test.select_atoms('resname CL')) == neg


@pytest.mark.parametrize('method, pos, neg', [
    ['ADD_NEUTRALIZE', 32, 25],
    ['SPLIT', 29, 22],
    ['SLTCAP', 28, 21],
    ['GENION', 36, 29]
])
def test_6olx(tmpdir, method, pos, neg):
    with tmpdir.as_cwd():
        genion(s=TPR_6OLX, structure=GRO_6OLX_NOIONS, o='test.gro',
               method=method)
        test = mda.Universe('test.gro')
        assert len(test.select_atoms('resname NA')) == pos
        assert len(test.select_atoms('resname CL')) == neg


def test_warnings_npnn_set(tmpdir):
    with tmpdir.as_cwd():
        with pytest.warns(UserWarning, match='nn passed to genion'):
            genion(s=TPR_6OLX, structure=GRO_6OLX_NOIONS, o='test.gro', nn=10)

        with pytest.warns(UserWarning, match='np passed to genion'):
            genion(s=TPR_6OLX, structure=GRO_6OLX_NOIONS, o='test.gro', np=0)


def test_pq_nq_typerror(tmpdir):
    with tmpdir.as_cwd():
        with pytest.raises(TypeError, match="Argument mismatch for pq"):
            genion(s=TPR_6OLX, structure=GRO_6OLX_NOIONS, o='test.gro', pq=10)

        with pytest.raises(TypeError, match="Argument mismatch for nq"):
            genion(s=TPR_6OLX, structure=GRO_6OLX_NOIONS, o='test.gro', nq=-5)
