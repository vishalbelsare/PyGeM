"""
PyGeM init
"""
# __all__ = [
#     'affine', 'filehandler', 'freeform', 'radial', 'openfhandler', 'elmerhandler',
#     'stlhandler', 'unvhandler', 'vtkhandler', 'nurbshandler', 'stephandler',
#     'igeshandler', 'utils', 'gui', 'khandler', 'idw'
# ]

def get_current_year():
    from datetime import datetime
    return datetime.now().year

__title__ = "pygem"
__author__ = "Marco Tezzele, Nicola Demo"
__copyright__ = "Copyright 2017-{}, PyGeM contributors".format(get_current_year())
__license__ = "MIT"
__version__ = "2.0.0"
__mail__ = 'marcotez@gmail.com, demo.nicola@gmail.com'
__maintainer__ = __author__
__status__ = "Stable"

#from .affine import *
from .deformation import Deformation
from .ffd import FFD
from .rbf import RBF
from .idw import IDW
from .rbf_factory import RBFFactory