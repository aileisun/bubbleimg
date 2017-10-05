# __init__.py
# ALS 2017/05/11

__all__ = ['obsobj', 'plainobj', 'sdss', 'hsc', 'objnaming', 'operator', 'imager']


import obsobj
import plainobj
import sdss
import hsc
import objnaming
import operator
import imager

# reload(obsobj)
# reload(plainobj)
# reload(sdss)
# reload(hsc)
# reload(objnaming)
# reload(operator)
# reload(imager)

from obsobj import obsObj
from operator import Operator
from imager import Imager

from sdss import sdssObj
from hsc import hscObj
