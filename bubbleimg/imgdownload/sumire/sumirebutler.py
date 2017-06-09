# sumirebutler.py
# ALS 2017/06/09

""" mimic the behavior of hscpipe butler to access data in the sumire cluster. It's initialized by rerun path, and offers methods datasetExists() and get(). Deal with gz decompression automatically """

import lsst.afw.image as afwImage


class Butler(object):

	def __init__(self):
		"""
		Butler
			sumirebutler that mimics lsst butler and runs in the sumire cluster. 

		Params
		------
		root (str)
			root to the rerun directory, e.g., '/array2/SSP/dr1/s16a/data/s16a_wide/'.

		dir_temp (str) ='/home/alsun/hscbbselect/temporary'
			directory to store 


		Public Methods
		--------------
		__init__(self, root)

		datasetExists(self, datasetType="deepCoadd_calexp", dataId={})

		get(self, datasetType="deepCoadd_calexp", dataId={})
		"""
		self.rootrerun_path = root
		self.dir_temp = dir_temp

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
		self.dir_temp = dir_temp
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
		fp = self._get_full_file_path(self, datasetType=datasetType, dataId=dataId)


		if os.path.isfile(fp):
			exp = afwImage.ExposureF(fp)
			return exp
		else:
			fp_gz = fp+'.gz'
			if os.path.isfile(fp_gz):

				with gzip.open(fp_gz, 'rb') as f_in, tempfile.NamedTemporaryFile(mode='w') as f_out:

					shutil.copyfileobj(f_in, f_out)
					exp = afwImage.ExposureF(f_out.name)
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
		fp = self._get_full_file_path(self, datasetType=datasetType, dataId=dataId)

		if os.path.isfile(fp):
			return True
		else:
			fp_gz = fp+'.gz'
			if os.path.isfile(fp_gz):
				return True
			else: 
				return False


# def get_decompressed_fits_filepath(fp):
# 	""" 
# 	take in a file path fp of a fits file that one wants to read. 
# 	if the file exist then return itself
# 	if not, check if its .gz exist, 
# 		if so, uncompress it and return a new path towards the uncompressed fits file
# 		otherwise, raise exception. 
# 	"""

# 	if os.path.isfile(fp):
# 		return fp
# 	else:
# 		fp_gz = fp+'.gz'
# 		if os.path.isfile(fp_gz):
# 			return True
# 		else: 
# 			raise Exception("[Butler] file does not exist")

# 	return fp_decomp


# def decompress_gz(fp_in):
# 	""" take in a file path fp_in and return a filepath towards an unzipped version of that file """

# 	if not os.path.isfile(fn_in):
# 		fn_in_gz = fn_in+'.gz'
# 		if os.path.isfile(fn_in_gz):


# 			ungz(fn_in_gz, fn_temp)
# 			fn_in_update = fn_temp
# 		else:
# 			raise Exception("[Butler] file does not exist")
# 	else:
# 		fn_in_update = fn_in

# 	return fp_out

# def ungz(fn_in, fn_out):
# 	""" decompress .gz file """
# 	with gzip.open(fn_in, 'rb') as f_in, open(fn_out, 'wb') as f_out:
# 		shutil.copyfileobj(f_in, f_out)

# import gzip
# import tempfile
# file = tempfile.NamedTemporaryFile()