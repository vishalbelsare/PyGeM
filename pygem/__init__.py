"""
PyGeM init
"""
# __all__ = [
#     'affine', 'filehandler', 'freeform', 'radial', 'openfhandler',
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

from .affine import *
from .freeform import FFD
from .radial import RBF
from .idw import IDW
from .filehandler import FileHandler
from .openfhandler import OpenFoamHandler
from .stlhandler import StlHandler
from .unvhandler import UnvHandler
from .vtkhandler import VtkHandler
from .khandler import KHandler
from .mdpahandler import MdpaHandler
from .params import *
