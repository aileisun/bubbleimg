# hscimgloader.py 
# ALS 2017/05/02

import numpy as np
import os
import requests
import astropy.units as u
from astropy.io import fits
import re

from ..loader import imgLoader
from ...filters import surveysetup
from ..get_credential import getCrendential
import hscurl
# from downloadbutler import multiButler
# import psf

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy


class hscimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		hscimgLoader, child of imgLoader

		download stamps from HSC DAS Quarry
		download psf by either: 
			(iaa) call sumire to infer psf from calexp and download psf from sumire
			(online) download calexp from server and infer psf locally

		on top of imgLoader init, set self.survey = 'hsc', 
		add attributes self.img_width_pix, self.img_height_pix
		do not load obj.sdss.xid by default unless to_make_obj_sdss= True


		Additional Params
		-----------------
		rerun = 's16a_wide' (str)
		release_version = 'dr1' (str)
		username (optional) (str): STARs account
		password (optional) (str): STARs account


		Public Methods
		--------------
		__init__(self, **kwargs)

		make_stamp(self, band, overwrite=False, **kwargs)

		make_stamps(self, overwrite=False, **kwargs)

		make_psf(self, band, overwrite=False, to_keep_calexp=False)

		make_psfs(self, overwrite=False, to_keep_calexp=False)


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


		Attributes
		----------
		(in addition to loader attributes)

		rerun = s16a_wide
		semester = s16a
		release_version = dr1
		survey = 'hsc'
		bands = ['g', 'r', 'i', 'z', 'y']
		username
		password
		status (bool)
			whether an hsc object is successfully identified

		"""
		super(hscimgLoader, self).__init__(**kwargs)

		# set data release parameters
		self.rerun = kwargs.pop('rerun', 's16a_wide')
		self.semester = self.rerun.split('_')[0]
		self.release_version = kwargs.pop('release_version', 'dr1')

		# set hsc object parameters
		self.status = super(self.__class__, self).add_obj_hsc(update=False, release_version=self.release_version, rerun=self.rerun)
		self.survey = 'hsc'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()

		# set connection parameters
		self.__username = kwargs.pop('username', '')
		self.__password = kwargs.pop('password', '')
		if self.__username == '' or self.__password == '':
			self.__username = getCrendential("HSC_SSP_CAS_USERNAME", cred_name = 'STARs username')
			self.__password = getCrendential("HSC_SSP_CAS_PASSWORD", cred_name = 'STARs password')


	def _get_fn_calexp(self, band):
		return 'calexp-{0}.fits'.format(band)


	def _get_filter_name(self, band):
		return "HSC-{band}".format(band=band.upper())


	def make_stamp(self, band, overwrite=False, **kwargs):
		"""
		make stamp image of the specified band of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. See _download_stamp() for specific implementation. 

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False
		**kwargs: to be passed to _download_stamp()
			e.g., imgtype='coadd', tract='', rerun='', see _download_stamp()
			if not specified then use self.rerun. 

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""
		return self._imgLoader__make_file_core(func_download_file=self._download_stamp, func_naming_file=self.get_fn_stamp, band=band, overwrite=overwrite, **kwargs)


	def make_stamps(self, overwrite=False, **kwargs):
		"""
		make stamps of all bands, see make_stamp()
		"""
		return self._imgLoader__make_files_core(func_download_file=self._download_stamp, func_naming_file=self.get_fn_stamp, overwrite=overwrite, **kwargs)


	def _download_stamp(self, band, imgtype='coadd', tract='', tokeepraw=False, n_trials=5):
		"""
		download hsc cutout img using HSC DAS Querry. Provides only ra, dec to DAS Querry and download the default coadd. always overwrite. 

		 convert it to stamp images. 

		ra, dec can be decimal degrees (12.345) or sexagesimal (1:23:35)

		for details see hsc query manual
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 

		Args
		--------
		band
		imgtype='coadd'
		tract=''
		tokeepraw = False (bool): 
			whether to keep the downloaded raw HSC image, which has four extensions. 
		n_trials=5
			how many times to retry requesting if there is requests errors such as connection error. 


		Return
		----------
		status: True if downloaded, False if download fails
		"""
		rerun = self.rerun

		# setting 
		fp_out = self.get_fp_stamp(band)

		semi_width_inarcsec = (self.img_width_arcsec.to(u.arcsec).value/2.)-0.1 # to get pix number right
		semi_height_inarcsec = (self.img_height_arcsec.to(u.arcsec).value/2.)-0.1
		sw = '%.5f'%semi_width_inarcsec+'asec'
		sh = '%.5f'%semi_height_inarcsec+'asec'

		# get url
		url = hscurl.get_hsc_cutout_url(self.ra, self.dec, band=band, rerun=rerun, tract=tract, imgtype=imgtype, sw=sw, sh=sh)

		# query, download, and convert to new unit
		# writing two files (if successful): raw img file and stamp img file. 
		rqst = self._retry_request(url, n_trials=n_trials)

		if rqst.status_code == 200:
			fp_raw = self._write_request_to_file(rqst)
			self._write_fits_unit_specified_in_nanomaggy(filein=fp_raw, fileout=fp_out)

			if not tokeepraw:
				os.remove(fp_raw)

			return True
		else:  
			print "[hscimgloader] image cannot be retrieved"
			return False


	def make_psf(self, band, overwrite=False, **kwargs):
		"""
		make psf image of the specified band of the object.  See _download_psf() for details. 

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False
		**kwargs: to be passed to _download_psf()
			e.g., imgtype='coadd'

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""
		return self._imgLoader__make_file_core(func_download_file=self._download_psf, func_naming_file=self.get_fn_psf, band=band, overwrite=overwrite, **kwargs)


	def make_psfs(self, overwrite=False, **kwargs):
		"""
		make psfs of all bands, see make_psf()
		"""
		return self._imgLoader__make_files_core(func_download_file=self._download_psf, func_naming_file=self.get_fn_psf, overwrite=overwrite, **kwargs)


	def _download_psf(self, band, imgtype='coadd', rerun='', tract='', patch_s='', n_trials=5):
		"""
		download hsc cutout img using HSC DAS Querry. Provides only ra, dec to DAS Querry and download the default coadd. always overwrite. If rerun not specified then use self.rerun. 

		for details see manual https://hscdata.mtk.nao.ac.jp/psf/4/manual.html#Bulk_mode
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 

		Args
		--------
		band
		imgtype='coadd'
		rerun=self.rerun
		tract=''
		patch_s=''
		n_trials=5
			how many times to retry requesting if there is requests errors such as connection error. 

		Return
		----------
		status: True if downloaded, False if download fails
		"""
		if rerun == '':
			rerun = self.rerun

		# setting 
		fp_out = self.get_fp_psf(band)

		# get url
		url = hscurl.get_hsc_psf_url(ra=self.ra, dec=self.dec, band=band, rerun=rerun, tract=tract, patch=patch_s, imgtype=imgtype)

		# download
		rqst = self._retry_request(url, n_trials=n_trials)

		if rqst.status_code == 200:
			self._write_request_to_file(rqst, fn=os.path.basename(fp_out))

			return True
		else:  
			print("[hscimgloader] psf cannot be retrieved")
			return False


	def _retry_request(self, url, n_trials=5):
		"""
		request url and retries for up to n_trials times if requests exceptions are raised, such as ConnectionErrors. Uses self.__username self.__password as authentication. 
		"""
		for _ in range(n_trials):
			try:
				rqst = requests.get(url, auth=(self.__username, self.__password))
				return rqst
				break
			except requests.exceptions.RequestException as e:
				print("[hscimgloader] retrying as error detected: "+str(e))


	def _write_request_to_file(self, rqst, fn=''):
		""" 
		write requested file under self.dir_obj with original filename unless filename specified

		Args
		--------
		rqst: request result
		fn ='' (str):
			the filename to be saved to. default: use original filename. 

		Return
		--------
		fp_out (string): the entire filepath to the file written
		"""
		d = rqst.headers['content-disposition']

		if fn == '':
			fn = re.findall("filename=(.+)", d)[0][1:-1]

		fp_out = self.dir_obj + fn

		with open(fp_out, 'wb') as out:
			for bits in rqst.iter_content():
				out.write(bits)
		return fp_out




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

