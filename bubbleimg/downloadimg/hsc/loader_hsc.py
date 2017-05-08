# loader.py 
# ALS 2017/05/02

import numpy as np
import os
import requests
import astropy.units as u
from astropy.io import fits
import re
import getpass


from ..loader import imgLoader
from ...filters import surveysetup

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy


class HSCimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		on top of imgLoader init, set self.survey = 'hsc', 
		add attributes self.img_width_pix, self.img_height_pix
		do not load obj.sdss.xid by default unless to_make_obj_sdss= True
		"""
		super(HSCimgLoader, self ).__init__(**kwargs)

		self.survey = 'hsc'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()

		if self.user == "" or self.password == "":
			self.user = raw_input("STARs username:")
			self.password = getpass.getpass("Password:")


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
		return self._imgLoader__make_stamp_core(self._download_stamp, band=band, overwrite=overwrite, **kwargs)


	def make_stamps(self, overwrite=False, **kwargs):
		return self._imgLoader__make_stamps_core(self._download_stamp, overwrite=overwrite, **kwargs)


	def _download_stamp(self, band='r', imgtype='coadd', tract='', rerun=''):
		"""
		download hsc cutout image

		ra, dec can be decimal degrees (12.345) or sexagesimal (1:23:35)

		for details see hsc query manual
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 

		Args
		--------
		band='r'
		imgtype='coadd'
		tract=''
		rerun=''

		Return
		----------
		status: True if downloaded, False if download fails

		!!! WARNING !!! 
		The original fits downloaded have four hdu, but only the second one is useful. To go around this problem. I save it and open it as astropy fits, keeping the image hdu and save it again. There are probably smarter ways of doing this without saving and reading. 
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
		rqst = requests.get(url, auth=(self.user, self.password))
		if rqst.status_code == 200:
			filepath_raw = self._write_request_to_file(rqst)
			self._write_fits_unit_converted_to_nanomaggy(filein=filepath_raw, fileout=filepath_out)

			return True
		else:  
			print "image cannot be retrieved"
			return False


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


	def _make_hsc_cutout_url(sefl, ra, dec, band='i', sw='5asec', sh='5asec', imgtype='coadd', tract='', rerun=''):
		"""
		see hsc query manual
		https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 
		"""

		url = 'https://hscdata.mtk.nao.ac.jp:4443/das_quarry/cgi-bin/quarryImage?ra={0}&dec={1}&sw={2}&sh={3}&type={4}&image=on&filter=HSC-{5}&tract={6}&rerun={7}'.format(ra, dec, sw, sh, imgtype, band.capitalize(), tract, rerun)

		return url



	def _write_fits_unit_converted_to_nanomaggy(self, filein, fileout):
		"""
		convert raw hsc image to stamp image
		take the second hdu hdu[1] as data

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
