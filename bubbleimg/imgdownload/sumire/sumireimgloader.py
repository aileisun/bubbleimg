# sumireimgloader.py 
# ALS 2017/06/09

import numpy as np
import os
# import requests
import astropy.units as u
from astropy.io import fits
# import re

from ..loader import imgLoader
from ...filters import surveysetup
from ..get_credential import getCrendential

import cutout
# from lsst.daf.persistence import Butler
from sumirebutler import Butler

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy


class sumireimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		sumireimgLoader, child of imgLoader

		Load hsc images within the iaa sumire cluster using lsst stack. Create calexp cutouts instead of plain stamp. 

		on top of imgLoader init, set self.survey = 'hsc', 
		add attributes self.img_width_pix, self.img_height_pix


		Additional Params
		-----------------
		rerun = 's16a_wide' (str)
		release_version = 'dr1' (str)

		"""
		super(sumireimgLoader, self).__init__(**kwargs)

		# set data release parameters
		self.rerun = kwargs.pop('rerun', 's16a_wide')
		self.semester = self.rerun.split('_')[0]
		self.release_version = kwargs.pop('release_version', 'dr1')
		self.rootrerun_path = '/array2/SSP/' + self._get_rerun_path()

		# set hsc object parameters
		self.status = super(self.__class__, self).add_obj_hsc(update=False, release_version=self.release_version, rerun=self.rerun)
		self.survey = 'hsc'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()


	def _get_fn_calexp(self, band):
		return 'calexp-{0}.fits'.format(band)


	def _get_filter_name(self, band):
		return "HSC-{band}".format(band=band.upper())


	def _get_rerun_path(self):
		"""
		e.g., "dr1/s16a/data/s16a_wide/"
		"""
		return '{dr}/{semester}/data/{rerun}/'.format(dr=self.release_version, semester=self.semester, rerun=self.rerun)

	def make_stamp(self, band, overwrite=False):
		"""
		make stamp image of the specified band of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. The stamps so creates are ExposureF files containing the mask, variance, psf header, etc.

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False
		# **kwargs: to be passed to _download_stamp()
		# 	e.g., imgtype='coadd', tract='', rerun='', see _download_stamp()

		Return
		----------
		status (bool)
		"""

		return self._imgLoader__make_file_core(func_download_file=self._download_stamp, func_naming_file=self.get_fn_stamp, band=band, overwrite=overwrite)


	def make_stamps(self, overwrite=False):
		"""
		make stamps of all bands, see make_stamp()
		"""
		return self._imgLoader__make_files_core(func_download_file=self._download_stamp, func_naming_file=self.get_fn_stamp, overwrite=overwrite)


	def _download_stamp(self, band, datasetType="deepCoadd_calexp"):
		"""
		make hsc calexp image cutouts and save to stamp-{band}.fits
		using hsc_xid.csv tract patch_s info to download, and ra, dec, img_width, etc to make cutouts. 
		always overwrites. 
		
		Args
		--------
		band (str)
		datasetType = "deepCoadd_calexp" (str)

		Return
		----------
		status: True if downloaded, False if download fails		
		"""
		fp = self.get_fp_stamp(band)

		b = Butler(self.rootrerun_path)

		dataId = dict(tract=self.obj.hsc.tract, patch=self.obj.hsc.patch_s, filter=self._get_filter_name(band))

		if b.datasetExists(datasetType, dataId):

			exp = b.get(datasetType, dataId)
			stamp = cutout.make_cutout(exp, self.ra, self.dec, npix=self.img_width_pix)
			stamp.writeFits(fp)
			return True
		else:
			print "[sumireimgloader] image cannot be retrieved"
			return False




	# taken from the old hscimgloader:
	# def make_psf(self, band, overwrite=False, to_keep_calexp=False):
	# 	"""
	# 	make psf files by downloading calexp file using selected method depending on self.environment and read psf info from it using psf.py

	# 	Params
	# 	------
	# 	band (str): e.g., 'r'
	# 	overwrite=False: whether to overwrite existing psf file
	# 	to_keep_calexp=False: whether to keep the calexp file or not after psf files are made
	# 	"""
	# 	if self.environment != 'iaa':
	# 		status = self._imgLoader__make_file_core(func_download_file=self._calexp_to_psf, func_naming_file=self.get_fn_psf, band=band, overwrite=overwrite)

	# 		if not to_keep_calexp:
	# 			fn = self.dir_obj+self._get_fn_calexp(band=band)
	# 			if os.path.isfile(fn):
	# 				os.remove(fn)
	# 	else: 
	# 		status = self._imgLoader__make_file_core(func_download_file=self._download_psf_at_iaa, func_naming_file=self.get_fn_psf, band=band, overwrite=overwrite)

	# 	return status


	# def make_psfs(self, overwrite=False, to_keep_calexp=False):

	# 	if self.environment!='iaa':
	# 		statuss = np.ndarray(5, dtype=bool)
	# 		for i, band in enumerate(self.bands): 
	# 			statuss[i] = self.make_psf(band=band, overwrite=overwrite, to_keep_calexp=to_keep_calexp)
	# 		return all(statuss)
	# 	else: 
	# 		missings = [(not os.path.isfile(self.dir_obj+self.get_fn_psf(band))) for band in self.bands]
	# 		isfilesmissing = np.any(missings)

	# 		if isfilesmissing or overwrite:
	# 			print "[hscimgloader] running _download_psfs_at_iaa()"
	# 			status = self._download_psfs_at_iaa()
	# 			return status
	# 		else: 
	# 			print "[hscimgloader] skip _download_psfs_at_iaa() as file exists"
	# 			return True


	# def _download_psf_at_iaa(self, band):
	# 	if self.environment!='iaa':
	# 		raise Exception("[hscimgLoader] _download_psf_at_iaa() can only be called when environment is iaa")

	# 	fn_out = self.dir_obj+self.get_fn_psf(band=band)

	# 	dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self._get_filter_name(band))

	# 	b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

	# 	with b:
	# 		status = b.download_psf(fn_out, ra=self.ra, dec=self.dec, **dataId)
	# 	return status


	# def _download_psfs_at_iaa(self):
	# 	if self.environment!='iaa':
	# 		raise Exception("[hscimgLoader] _download_psf_at_iaa() can only be called when environment is iaa")

	# 	b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

	# 	statuss = np.ndarray(len(self.bands), dtype=bool)
	# 	with b:
	# 		for i, band in enumerate(self.bands): 
	# 			fn_out = self.dir_obj+self.get_fn_psf(band=band)

	# 			dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self._get_filter_name(band))

	# 			statuss[i] = b.download_psf(fn_out, ra=self.ra, dec=self.dec, **dataId)
	# 	return all(statuss)


	# def _calexp_to_psf(self, band):
	# 	fn_in = self.dir_obj+self._get_fn_calexp(band=band)
	# 	fn_out = self.dir_obj+self.get_fn_psf(band=band)

	# 	status_calexp = self._make_calexp(band=band, overwrite=False)

	# 	if status_calexp:
	# 		psf.exposureF_to_psf(fn_in, fn_out, self.ra, self.dec)

	# 		status = os.path.isfile(fn_out)
	# 		return status
	# 	else: False


	# def _make_calexp(self, band, overwrite=False):
	# 	"""
	# 	make calexp image of the specified band of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

	# 	Params
	# 	----------
	# 	band (string)
	# 	overwrite (boolean) = False

	# 	Return
	# 	----------
	# 	status: True if downloaded or skipped, False if download fails
	# 	"""
	# 	status = self._imgLoader__make_file_core(func_download_file=self._download_calexp, func_naming_file=self._get_fn_calexp, band=band, overwrite=overwrite)

	# 	return status


	# def _make_calexps(self, overwrite=False):
	# 	"""
	# 	make calexps of all bands, see _make_calexp()
	# 	"""
	# 	status = self._imgLoader__make_files_core(func_download_file=self._download_calexp, func_naming_file=self._get_fn_calexp, overwrite=overwrite)

	# 	return status


	# def _download_calexp(self, band):

	# 	fn = self.dir_obj+self._get_fn_calexp(band=band)

	# 	dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self._get_filter_name(band))

	# 	b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

	# 	with b:
	# 		status = b.download_file(fn, **dataId)
	# 	return status
