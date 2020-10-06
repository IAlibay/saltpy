"""
saltpy
A python toolset for calculating and adding salt concentration to solvated atomistic systems
"""

# Add imports here
# from .calculate import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
