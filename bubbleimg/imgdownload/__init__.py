# __init__.py
# ALS 2017/05/02

__all__ = ['loader']

from . import loader, sdss, hsc

from .loader import imgLoader
from .sdss.sdssimgloader import sdssimgLoader
from .hsc.hscimgloader import hscimgLoader