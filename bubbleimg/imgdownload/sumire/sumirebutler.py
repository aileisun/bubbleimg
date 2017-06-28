# sumirebutler.py
# ALS 2017/06/09

""" mimic the behavior of hscpipe butler to access data in the sumire cluster. It's initialized by rerun path, and offers methods datasetExists() and get(). Deal with gz decompression automatically """

import os
import gzip
import tempfile
import shutil

import lsst.afw.image as afwImage

from astropy.io import fits

class Butler(object):

	def __init__(self, root):
		"""
		Butler
			sumirebutler that mimics lsst butler and runs in the sumire cluster. 

		Params
		------
		root (str)
			root to the rerun directory, e.g., '/array2/SSP/dr1/s16a/data/s16a_wide/'.

		Public Methods
		--------------
		__init__(self, root)

		datasetExists(self, datasetType="deepCoadd_calexp", dataId={})

		get(self, datasetType="deepCoadd_calexp", dataId={})
		"""
		self.rootrerun_path = root


	def _get_tail_path(self, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'):
		"""
		Return
		------
		tail_path (str) :
			e.g., "deepCoadd/HSC-R/9564/7,3/calexp-HSC-R-9564-7,3.fits"
		"""
		tail_path = '{coadd}/{filter}/{tract}/{patch_s}/{filetype}-{filter}-{tract}-{patch_s}.fits'.format(coadd=coadd, filetype=filetype, tract=str(tract), patch_s=str(patch_s), filter=filter)
		return tail_path


	def _get_full_file_path(self, datasetType="deepCoadd_calexp", dataId={}):
		"""
		Params
		------
		datasetType = "deepCoadd_calexp"
		dataId = {}:
			list of dataid
			e.g., dict(tract=9564, patch="7,3", filter="HSC-R")

		Return
		------
		fp (str)
		"""
		coadd, filetype = datasetType.split('_')
		tract, patch_s, filter = dataId['tract'], dataId['patch'], dataId['filter']
		fp = self.rootrerun_path + self._get_tail_path(tract, patch_s, filter, coadd=coadd, filetype=filetype)
		return fp


	def get(self, datasetType="deepCoadd_calexp", dataId={}):
		"""
		Params
		------
		datasetType = "deepCoadd_calexp"
		dataId = {}:
			list of dataid
			e.g., dict(tract=9564, patch="7,3", filter="HSC-R")

		Return
		------
		exp
		"""
		fp = self._get_full_file_path(datasetType=datasetType, dataId=dataId)


		if os.path.isfile(fp):
			exp = afwImage.ExposureF(fp)
			return exp
		else:
			fp_gz = fp+'.gz'
			if os.path.isfile(fp_gz):

				with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f_temp:
					with gzip.open(fp_gz, 'rb') as f_in:
						shutil.copyfileobj(f_in, f_temp)

					f_temp.close()

					exp = afwImage.ExposureF(f_temp.name)
					
					if os.path.isfile(f_temp.name):
						os.remove(f_temp.name)

					return exp

			else: 
				return Exception("[Butler] file does not exist")


	def datasetExists(self, datasetType="deepCoadd_calexp", dataId={}):
		"""
		Params
		------
		datasetType = "deepCoadd_calexp"
		dataId = {}:
			list of dataid
			e.g., dict(tract=9564, patch="7,3", filter="HSC-R")

		Return
		------
		result (bool)
		"""
		fp = self._get_full_file_path(datasetType=datasetType, dataId=dataId)

		if os.path.isfile(fp):
			return True
		else:
			fp_gz = fp+'.gz'
			if os.path.isfile(fp_gz):
				return True
			else: 
				return False
