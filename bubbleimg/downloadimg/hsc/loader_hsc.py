# loader.py 
# ALS 2017/05/02

import numpy as np
import os
import requests
import astropy.units as u
from astropy.io import fits
import re

from ..loader import imgLoader
from ...filters import surveysetup
from multibutler import multiButler
from get_credential import getCrendential
import psf

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy


class hscimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		hscimgLoader, child of imgLoader

		on top of imgLoader init, set self.survey = 'hsc', 
		add attributes self.img_width_pix, self.img_height_pix
		do not load obj.sdss.xid by default unless to_make_obj_sdss= True


		Additional Params
		-----------------
		environment = 'iaa' (str)
			: 'iaa', 'sumire', or 'online', chooses how large files are downloaded, see multibutler.py
		rerun = 's16a_wide' (str)
		release_version = 'dr1' (str)
		username (optional) (str): STARs account
		password (optional) (str): STARs account


		Instruction for stars username and password
		-------------------------------------------
		1) as arguments 
			hscimgLoader(..., username=username, password=password)

		2) as environmental variable 
		  $ export HSC_SSP_CAS_USERNAME
		  $ read -s HSC_SSP_CAS_USERNAME
		  $ export HSC_SSP_CAS_PASSWORD
		  $ read -s HSC_SSP_CAS_PASSWORD

		3) enter from terminal

		"""
		super(hscimgLoader, self).__init__(**kwargs)

		self.environment = kwargs.pop('environment', 'iaa')
		self.rerun = kwargs.pop('rerun', 's16a_wide')
		self.semester = self.rerun.split('_')[0]
		self.release_version = kwargs.pop('release_version', 'dr1')

		self.hsc_status = super(self.__class__, self).add_obj_hsc(update=False, release_version=self.release_version, rerun=self.rerun)

		self.survey = 'hsc'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()

		# get stars username ans password
		self.__username = kwargs.pop('username', '')
		self.__password = kwargs.pop('password', '')
		if self.__username == '' or self.__password == '':
			self.__username = getCrendential("HSC_SSP_CAS_USERNAME", cred_name = 'STARs username')
			self.__password = getCrendential("HSC_SSP_CAS_PASSWORD", cred_name = 'STARs password')


	def make_stamp(self, band='r', overwrite=False, **kwargs):
		"""
		make stamp image of the specified band of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False
		**kwargs: to be passed to _download_stamp()
			e.g., imgtype='coadd', tract='', rerun='', see _download_stamp()

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""

		# return self._imgLoader__make_stamp_core(func_download_stamp=self._download_stamp, band=band, overwrite=overwrite, **kwargs)

		return self._imgLoader__make_file_core(func_download_file=self._download_stamp, func_naming_file=self.get_stamp_filename, band=band, overwrite=overwrite, **kwargs)


	def make_stamps(self, overwrite=False, **kwargs):
		"""
		make stamps of all bands, see make_stamp()
		"""
		return self._imgLoader__make_files_core(func_download_file=self._download_stamp, func_naming_file=self.get_stamp_filename, overwrite=overwrite, **kwargs)


	def _download_stamp(self, band='r', imgtype='coadd', tract='', rerun='', tokeepraw=False):
		"""
		download hsc cutout image and convert it to stamp images. 

		ra, dec can be decimal degrees (12.345) or sexagesimal (1:23:35)

		for details see hsc query manual
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 

		Args
		--------
		band='r'
		imgtype='coadd'
		tract=''
		rerun=''
		tokeepraw=False (bool): whether to keep the downloaded raw HSC image, which has four extensions. 

		Return
		----------
		status: True if downloaded, False if download fails

		"""

		# setting 
		filepath_out = self.get_stamp_filepath(band)

		
		semi_width_inarcsec = (self.img_width_arcsec.to(u.arcsec).value/2.)-0.1 # to get pix number right
		semi_height_inarcsec = (self.img_height_arcsec.to(u.arcsec).value/2.)-0.1
		sw = '%.5f'%semi_width_inarcsec+'asec'
		sh = '%.5f'%semi_height_inarcsec+'asec'

		# get url
		url = self._make_hsc_cutout_url(self.ra, self.dec, band=band, sw=sw, sh=sh, imgtype=imgtype, tract=tract, rerun=rerun)

		# query, download, and convert to new unit
		# writing two files (if successful): raw img file and stamp img file. 
		rqst = requests.get(url, auth=(self.__username, self.__password))
		if rqst.status_code == 200:
			filepath_raw = self._write_request_to_file(rqst)
			self._write_fits_unit_specified_in_nanomaggy(filein=filepath_raw, fileout=filepath_out)

			if not tokeepraw:
				os.remove(filepath_raw)

			return True
		else:  
			print "image cannot be retrieved"
			return False


	def make_psf(self, band, overwrite=False, to_keep_calexp=False):
		"""
		make psf files by downloading calexp file using selected method depending on self.environment and read psf info from it using psf.py

		Params
		------
		band (str): e.g., 'r'
		overwrite=False: whether to overwrite existing psf file
		to_keep_calexp=False: whether to keep the calexp file or not after psf files are made
		"""
		if self.environment!='iaa':
			status = self._imgLoader__make_file_core(func_download_file=self._calexp_to_psf, func_naming_file=self.get_psf_filename, band=band, overwrite=overwrite)

			if not to_keep_calexp:
				fn = self.dir_obj+self.get_calexp_filename(band=band)
				if os.path.isfile(fn):
					os.remove(fn)
		else: 
			status = self._imgLoader__make_file_core(func_download_file=self._download_psf_at_iaa, func_naming_file=self.get_psf_filename, band=band, overwrite=overwrite)

		return status


	def make_psfs(self, overwrite=False, to_keep_calexp=False):

		if self.environment!='iaa':
			statuss = np.ndarray(5, dtype=bool)
			for i, band in enumerate(self.bands): 
				statuss[i] = self.make_psf(band=band, overwrite=overwrite, to_keep_calexp=to_keep_calexp)
			return all(statuss)
		else: 
			missings = [(not os.path.isfile(self.dir_obj+self.get_psf_filename(band))) for band in self.bands]
			isfilesmissing = np.any(missings)

			if isfilesmissing:
				return self._download_psfs_at_iaa()
			else: 
				print "skip _download_psfs_at_iaa() as file exists"
				return True


	def _download_psf_at_iaa(self, band):
		if self.environment!='iaa':
			raise Exception("[hscimgLoader] _download_psf_at_iaa() can only be called when environment is iaa")

		fn_out = self.dir_obj+self.get_psf_filename(band=band)

		dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self.get_filter_name(band))

		b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

		with b:
			status = b.download_psf(fn_out, ra=self.ra, dec=self.dec, **dataId)
		return status


	def _download_psfs_at_iaa(self):
		if self.environment!='iaa':
			raise Exception("[hscimgLoader] _download_psf_at_iaa() can only be called when environment is iaa")

		b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

		statuss = np.ndarray(len(self.bands), dtype=bool)
		with b:
			for i, band in enumerate(self.bands): 
				fn_out = self.dir_obj+self.get_psf_filename(band=band)

				dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self.get_filter_name(band))

				statuss[i] = b.download_psf(fn_out, ra=self.ra, dec=self.dec, **dataId)
		return all(statuss)


	def _calexp_to_psf(self, band):
		fn_in = self.dir_obj+self.get_calexp_filename(band=band)
		fn_out = self.dir_obj+self.get_psf_filename(band=band)

		status_calexp = self.make_calexp(band=band, overwrite=False)

		if status_calexp:
			psf.exposureF_to_psf(fn_in, fn_out, self.ra, self.dec)

			status = os.path.isfile(fn_out)
			return status
		else: False


	def make_calexp(self, band='r', overwrite=False):
		"""
		make calexp image of the specified band of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""
		return self._imgLoader__make_file_core(func_download_file=self._download_calexp, func_naming_file=self.get_calexp_filename, band=band, overwrite=overwrite)


	def make_calexps(self, overwrite=False):
		"""
		make calexps of all bands, see make_calexp()
		"""
		return self._imgLoader__make_files_core(func_download_file=self._download_calexp, func_naming_file=self.get_calexp_filename, overwrite=overwrite)


	def _download_calexp(self, band='r'):

		fn = self.dir_obj+self.get_calexp_filename(band=band)

		dataId = dict(tract=self.obj.hsc.tract, patch_s=self.obj.hsc.patch_s, filter=self.get_filter_name(band))

		b = multiButler(environment=self.environment, release_version=self.release_version, semester=self.semester, rerun=self.rerun).butler

		with b:
			status = b.download_file(fn, **dataId)
		return status




	def _write_request_to_file(self, rqst):
		""" 
		write requested file under self.dir_obj with original filename

		Args
		--------
		rqst: request result

		Return
		--------
		filepath_raw (string): the entire filepath to the file written
		"""
		d = rqst.headers['content-disposition']
		filename_raw = re.findall("filename=(.+)", d)[0][1:-1]
		filepath_raw = self.dir_obj + filename_raw

		with open(filepath_raw, 'wb') as out:
			for bits in rqst.iter_content():
				out.write(bits)
		return filepath_raw


	def _make_hsc_cutout_url(sefl, ra, dec, band='i', sw='5asec', sh='5asec', imgtype='coadd', tract='', rerun='', mask='on', variance='on'):
		"""
		see hsc query manual
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 
		"""

		url = 'https://hscdata.mtk.nao.ac.jp:4443/das_quarry/cgi-bin/quarryImage?ra={0}&dec={1}&sw={2}&sh={3}&type={4}&image=on&mask={5}&variance={6}&filter=HSC-{7}&tract={8}&rerun={9}'.format(ra, dec, sw, sh, imgtype, mask, variance, band.capitalize(), tract, rerun)

		return url


	def _write_fits_unit_converted_to_nanomaggy(self, filein, fileout):
		"""
		!!!!!!! WARNING !!!!!!!! this funciton is not used currently

		Convert raw hsc image to an image with unit nanomaggy, changing the data value. 
		take only the second hdu hdu[1] as data in output

		read in fits file filein with no bunit but FLUXMAG0 and convert to one fits file with unit nanomaggy, and write to fileout. 

		Notes on Unit conversion
		-----------
		HSC fluxmag0 is set such that a pix value of 1 has a magAB of 27 so:

			fluxmag0 = header_combine['FLUXMAG0']
			# 63095734448.0194

			pixunit = 10.**-19.44 / fluxmag0 * (u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
			# u.Quantity('5.754399373371546e-31 erg / cm2 / Hz / s')

			nanomaggy_per_raw_unit = float((u.nanomaggy/pixunit).decompose()) 
			# 63.099548091890085

		But this function should work even with other fluxmag 0, as we set 		
			nanomaggy_per_raw_unit = fluxmag0 * 10**-9

		"""
		hdu = fits.open(filein)
		header_combine = hdu[1].header+hdu[0].header

		# sanity check
		if header_combine['FLUXMAG0'] != 63095734448.0194: 
			raise ValueError("HSC FLUXMAG0 different from usual. Although shouldnt be a problem")
		if 'BUNIT' in header_combine: 
			raise ValueError("Input fits file should not have BUNIT")

		nanomaggy_per_raw_unit = header_combine['FLUXMAG0']*10.**-9
		data_nanomaggy = hdu[1].data/nanomaggy_per_raw_unit

		header_combine.set(keyword='BUNIT', value='nanomaggy', comment="1 nanomaggy = 3.631e-6 Jy")
		header_combine['COMMENT'] = "Unit converted to nanomaggy by ALS"
		header_combine.remove(keyword='FLUXMAG0')

		hdu_abbrv = fits.PrimaryHDU(data_nanomaggy, header=header_combine)
		hdu_abbrv.writeto(fileout, overwrite=True)


	def _write_fits_unit_specified_in_nanomaggy(self, filein, fileout):
		"""
		Convert a raw hsc image to an image with unit nanomaggy, the data values unchanged. 
		Take only the second hdu hdu[1] as data in output. 

		read in fits file filein with no bunit but FLUXMAG0 and convert to one fits file with unit nanomaggy, and write to fileout. 

		Notes on Unit conversion
		-----------
		HSC fluxmag0 is set such that a pix value of 1 has a magAB of 27 so:

			fluxmag0 = header_combine['FLUXMAG0']
			# 63095734448.0194

			pixunit = 10.**-19.44 / fluxmag0 * (u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
			# u.Quantity('5.754399373371546e-31 erg / cm2 / Hz / s')

			nanomaggy_per_raw_unit = float((u.nanomaggy/pixunit).decompose()) 
			# 63.099548091890085

			raw_unit_per_nanomaggy = 1/nanomaggy_per_raw_unit
			# 0.015847974038478506

		But this function should work even with other fluxmag 0, as we set 		
			nanomaggy_per_raw_unit = fluxmag0 * 10**-9

		"""
		hdu = fits.open(filein)
		header_combine = hdu[1].header+hdu[0].header

		# sanity check
		if header_combine['FLUXMAG0'] != 63095734448.0194: 
			raise ValueError("HSC FLUXMAG0 different from assumed")
		
		if 'BUNIT' in header_combine: 
			raise ValueError("Input fits file should not have BUNIT")
		
		bunit = '1.58479740e-02 nanomaggy'
		header_combine.set(keyword='BUNIT', value=bunit, comment="1 nanomaggy = 3.631e-6 Jy")
		header_combine['COMMENT'] = "Unit specified in nanomaggy by ALS"

		data = hdu[1].data
		hdu_abbrv = fits.PrimaryHDU(data, header=header_combine)
		hdu_abbrv.writeto(fileout, overwrite=True)


	def get_calexp_filename(self, band):
		return 'calexp-{0}.fits'.format(band)


	def get_filter_name(self, band):
		return "HSC-{band}".format(band=band.upper())