# obsobj.py
# ALS 2017/05/11

"""
define class obsObj 
"""

from plainobj import plainObj
from sdss.sdssobj import SDSSObj
from hsc.hscobj import HSCObj

class obsObj(plainObj):
	def __init__(self, **kwargs):
		"""
		obsObj

		Params
		----------
		ra (float)
		dec (float)
		/either
			dir_obj (string)
		/or 
			dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'
				

		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		"""

		super(self.__class__, self).__init__(**kwargs)


	def add_sdss(self, **kwargs):

		self.sdss = SDSSObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj, **kwargs)

		return self.sdss.status


	def add_hsc(self, **kwargs):
		self.hsc = HSCObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj, **kwargs)

		return self.hsc.status
