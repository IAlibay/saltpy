"""
Basic import test for the saltpy package.

Author: Irfan Alibay
Date: 2020
"""

# Import package, test suite, and other packages as needed
import saltpy
import pytest
import sys

def test_saltpy_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "saltpy" in sys.modules
