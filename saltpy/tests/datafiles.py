"""
Datafiles for tests
"""

__all__ = [
        "TPR_6OLX", "GRO_6OLX_NOIONS", "GRO_6OLX_IONS",
        "TPR_2RBN", "GRO_2RBN_NOIONS", "GRO_2RBN_IONS"
]

from pkg_resources import resource_filename

TPR_2RBN = resource_filename(__name__, '../data/2RBN/2RBN.cubic.genion.tpr')
GRO_2RBN_NOIONS = resource_filename(__name__,
                           '../data/2RBN/2RBN.cubic.gro')
GRO_2RBN_IONS = resource_filename(__name__,
                            '../data/2RBN/2RBN.cubic.150mM.gro')

TPR_6OLX = resource_filename(__name__,
                            '../data/6OLX/6OLX.dodeca.genion.tpr')
GRO_6OLX_NOIONS = resource_filename(__name__,
                            '../data/6OLX/6OLX.dodeca.gro')
GRO_6OLX_IONS = resource_filename(__name__,
                            '../data/6OLX/6OLX.dodeca.150mM.gro')
