"""
Basic import test for the saltpy package.

Author: Irfan Alibay
Date: 2020
"""

import sys


def test_saltpy_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "saltpy" in sys.modules
