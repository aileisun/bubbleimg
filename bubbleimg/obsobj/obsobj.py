# obsobj.py
# ALS 2017/05/11

from plainobj import plainObj
from sdss.sdssobj import sdssObj
from hsc.hscobj import hscObj


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


	def add_sdss(self, overwrite=False, **kwargs):

		self.sdss = sdssObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj, obj_naming_sys=self.obj_naming_sys, overwrite=overwrite, **kwargs)

		return self.sdss.status


	def add_hsc(self, overwrite=False, **kwargs):
		self.hsc = hscObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj, obj_naming_sys=self.obj_naming_sys, overwrite=overwrite, **kwargs)

		return self.hsc.status
