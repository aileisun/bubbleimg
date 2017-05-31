# hscbatch.py
# ALS 2017/05/29

import numpy as np

from ..batch import Batch
from ... import downloadimg
from ... import blobmaps

class hscBatch(Batch):
	def __init__(self, **kwargs): 
		kwargs['survey'] = 'hsc'
		super(self.__class__, self).__init__(**kwargs)


	def build(self, func_build=None, overwrite=False, **kwargs):
		"""
		build the batch

		Params
		------
		func_build=self._func_build (funcion):
			a funciton that takes (ra, dec, dir_parent, overwrite, **kwargs) as param and returns status (bool)
		overwrite=False (bool)
		**kwargs:
			 to be entered into func_build() in the kwargs part, e.g., 'environment'='iaa'. 

		Return 
		------
		status
		"""

		if func_build is None:
			func_build = self._func_build

		status = super(self.__class__, self)._batch__build_core(func_build, overwrite=overwrite, **kwargs)
		return status


	def _func_build(self, ra, dec, dir_parent, overwrite=False, **kwargs):
		"""
		Params
		------
		ra
		dec
		obj_name
		overwrite=False

		**kwargs:
			environment='iaa'

		Return
		------
		status
		"""

		# setting
		environment = kwargs.pop('environment', 'iaa')
		humvi_bands = 'riz'

		# running
		L = downloadimg.hscimgLoader(ra=ra, dec=dec, dir_parent=dir_parent, environment=environment)

		statuss = [
					L.hsc_status, 
					L.add_obj_sdss(), 
					L.make_stamps(overwrite=overwrite), 
					L.make_psfs(overwrite=overwrite), 
					L.obj.sdss.make_spec(overwrite=overwrite)
					]

		blobmaps.imagedisp_util.objw_HumVIgriimages(L, bands=humvi_bands, update=overwrite)

		return all(statuss)