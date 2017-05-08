# __init__.py
# ALS 2017/05/02

__all__ = ['loader']

import loader, sdss, hsc

reload(loader)
reload(sdss)
reload(hsc)

from loader import imgLoader
from sdss.loader_sdss import SDSSimgLoader
from hsc.loader_hsc import HSCimgLoader