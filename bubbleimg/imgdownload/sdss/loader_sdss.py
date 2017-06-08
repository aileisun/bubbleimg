# loader.py 
# ALS 2017/05/02

import numpy as np
import os
from astroquery.sdss import SDSS
from astropy.io import fits

from ..loader import imgLoader
from ...filters import surveysetup
import stamp
import psf

class sdssimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		on top of imgLoader init, set self.survey = 'sdss'
		force to_make_obj_sdss=True such that self.obj.sdss.xid is always loaded 
		add attributes self.img_width_pix, self.img_height_pix
		"""

		super(sdssimgLoader, self).__init__(**kwargs)
		self.status = super(self.__class__, self).add_obj_sdss(update=False)


		self.survey = 'sdss'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()

		# sanity check
		if self.status:
			if round(self.obj.sdss.xid['ra'][0], 4) != round(self.ra, 4):
				raise ValueError("xid created inconsistent with init ra")


	def make_stamps(self, overwrite=False, band_rf='r', tokeepframe=False):
		"""
		make stamp images of all the bands of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		overwrite (boolean) = False

		band_rf='r': reference band
		tokeepframe=True: if True to keep frame images, if False delete frame images

		Return
		------
		status: True if all stamp-(band).fits exist
		"""
		# make frame images
		for band in self.bands:
			self._make_frame(band=band, overwrite=overwrite)

		# make stamp images from frame
		isstampfiles = np.all([os.path.isfile(self.get_fp_stamp(b)) for b in self.bands])

		if (not isstampfiles) or overwrite:
			stamp.write_alignedstampImages(obj=self.obj, bands=self.bands, band_rf=band_rf, xwidth=self.img_width_pix, ywidth=self.img_height_pix, clipnegative=False, overwrite=True)

		if not tokeepframe:
			for band in self.bands:
				os.remove(self._get_frame_filepath(band))

		isstampfiles = np.all([os.path.isfile(self.get_fp_stamp(b)) for b in self.bands])
		return isstampfiles


	def get_stamp(self, band='r'):
		"""
		read and return stamp-(band).fits, and calls make_stamps if file does not exist. If all fail, returns false. 
		"""

		fn = self.get_fp_stamp(band=band)

		if os.path.isfile(fn):
			return fits.getdata(fn)
		else: 
			status = self.make_stamps()
			if status: 
				return fits.getdata(fn)
			else: 
				return False


	def make_psfs(self, overwrite=False, tokeepfield=False):
		"""
		make psfs of all the bands of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		overwrite (boolean) = False
		tokeepfield = True: if True to keep psField.fit, if False delete psField.fit

		Return
		------
		status: True if all psf-(band).fits exist
		"""
		# make psField.fit
		filename_psfield = 'psField.fits'
		if not os.path.isfile(self.dir_obj+filename_psfield) or overwrite:
			psf.download_psField(xid=self.obj.sdss.xid, dir_out=self.dir_obj, filename=filename_psfield)

	
		# make psf-(band).fit from psField.fit
		ispsffiles = np.all([os.path.isfile(self.get_fp_psf(b)) for b in self.bands])

		if (not ispsffiles) or overwrite:
			psf.psField_to_psfs(dir_obj=self.dir_obj, photoobj=self.obj.sdss.photoobj, bands=['u', 'g', 'r', 'i', 'z'])

		if not tokeepfield:
			os.remove(self.dir_obj+filename_psfield)

		ispsffiles = np.all([os.path.isfile(self.get_fp_psf(b)) for b in self.bands])
		return ispsffiles
		

	def get_psf(self, band='r'):
		"""
		read and return psf-(band).fits, and calls make_psfs if file does not exist. If all fail, returns false. 
		"""

		fn = self.get_fp_psf(band=band)

		if os.path.isfile(fn):
			return fits.getdata(fn)
		else: 
			status = self.make_psfs()
			if status: 
				return fits.getdata(fn)
			else: 
				return False


	def _make_frame(self, band='r', overwrite=True):
		""" 
		download frame image if it does not exist

		Params
		-------
		band (string) = 'r'

		Return
		-------
		status (bool): True if successful, False if not
		"""
		file = self._get_frame_filepath(band)

		if (not os.path.isfile(file)) or overwrite:
			print "[sdssimgLoader] downloading frame image band {0}".format(band)
			im = SDSS.get_images(matches=self.obj.sdss.xid, band=band)
			im[0].writeto(file, overwrite=True)
		else: 
			print "[sdssimgLoader] skip downloading frame image band {0}".format(band)

		return os.path.isfile(file)


	def _get_frame_filepath(self, band):
		return self.dir_obj+'frame-{0}.fits'.format(band)


