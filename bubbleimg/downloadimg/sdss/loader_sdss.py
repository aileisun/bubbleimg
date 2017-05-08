# loader.py 
# ALS 2017/05/02

import numpy as np
import os
from astroquery.sdss import SDSS

from ..loader import imgLoader
from ...filters import surveysetup
import alignstamp

class SDSSimgLoader(imgLoader):

	def __init__(self, **kwargs):
		""" 
		on top of imgLoader init, set self.survey = 'sdss'
		force to_make_obj_sdss=True such that self.obj.sdss.xid is always loaded 
		add attributes self.img_width_pix, self.img_height_pix
		"""
		kwargs['to_make_obj_sdss'] = True
		super(SDSSimgLoader, self ).__init__(**kwargs)

		# make sure there is self.obj.sdss.xid
		self.survey = 'sdss'
		self.bands = surveysetup.surveybands[self.survey]
		self.pixsize = surveysetup.pixsize[self.survey]
		self._add_attr_img_width_pix_arcsec()

		# sanity check
		if round(self.obj.sdss.xid['ra'][0], 4) != round(self.ra, 4):
			raise ValueError("xid created inconsistent with init ra")


	def make_stamps(self, overwrite=False, band_rf='r', tokeepframe=False):
		"""
		make stamp images of all the bands of the object. takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		band (string) = 'r'
		overwrite (boolean) = False

		band_rf='r': reference band
		tokeepframe=True: if True to keep frame images, if False delete frame images
		"""
		for band in self.bands:
			self._make_frame(band=band)

		isstampfiles = np.all([os.path.isfile(self.get_stamp_filepath(b)) for b in self.bands])

		if (not isstampfiles) or overwrite:
			alignstamp.getalignedstampImages(obj=self.obj, bands=self.bands, band_rf=band_rf, xwidth=self.img_width_pix.value, ywidth=self.img_height_pix.value, savefits=True, clipnegative=False)

		if not tokeepframe:
			for band in self.bands:
				os.remove(self.get_frame_filepath(band))

		isstampfiles = np.all([os.path.isfile(self.get_stamp_filepath(b)) for b in self.bands])
		return isstampfiles


	def _make_frame(self, band='r'):
		""" 
		download frame image if it does not exist
		WARNING: frame image is never overwritten

		Params
		-------
		band (string) = 'r'

		Return
		-------
		status (bool): True if successful, False if not
		"""
		file = self.get_frame_filepath(band)

		if not os.path.isfile(file):
			print "[SDSSimgLoader] downloading frame image band {0}".format(band)
			im = SDSS.get_images(matches=self.obj.sdss.xid, band=band)
			im[0].writeto(file)
		else: 
			print "[SDSSimgLoader] skip downloading frame image band {0}".format(band)

		return os.path.isfile(file)


	def get_frame_filepath(self, band):
		return self.dir_obj+'frame-{0}.fits'.format(band)


