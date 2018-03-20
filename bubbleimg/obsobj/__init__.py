# __init__.py
# ALS 2017/05/11

__all__ = ['obsobj', 'plainobj', 'sdss', 'hsc', 'objnaming', 'operator', 'imager']


from . import obsobj
from . import plainobj
from . import sdss
from . import hsc
from . import objnaming
from . import operator
from . import imager

from .obsobj import obsObj
from .operator import Operator
from .imager import Imager

from .sdss import sdssObj
from .hsc import hscObj
