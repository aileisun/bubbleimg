# __init__.py
# ALS 2017/05/11

__all__ = ['sdssobj']

from . import sdssobj
import imp

imp.reload(sdssobj)

from .sdssobj import sdssObj
