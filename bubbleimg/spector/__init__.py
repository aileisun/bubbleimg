# __init__.py

__all__ = ['spector']

from . import spector
import imp
# import linefrac

imp.reload(spector)
# reload(linefrac)

from .spector import Spector
