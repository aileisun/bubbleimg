# decomposer.py 
# ALS 2017/06/01

import os

from ..obsobj import Imager
from .. import spector
from .. import imgdownload
from ..filters import surveysetup


class Decomposer(Imager):

	def __init__(self, **kwargs):
		"""
		Decomposer, an imager operator to do image decomposition

		Imager Params
		-------------
		/either
			obj (object of class obsobj): with attributes ra, dec, dir_obj
		/or  
			ra (float)
			dec (float)
			/either
				dir_obj (string)
			/or 
				dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'

		survey (str): 
			survey of the photometric system
			if not provided, use self.obj.survey. Raise exception if self.obj.survey does not exist. 
		z (float):  
			redshift, if not provided, use self.obj.z or self.obj.sdss.z. It does not automatically query sdss to get z. If nothing is pecified then set to -1. 

		center_mode='n/2' (str):
			how is image center defined in case of even n, 'n/2' or 'n/2-1'. Should be set to n/2-1 if the image is downloaded from HSC quarry. 


		Imager Attributes
		-----------------
		obj (instance of objObj)
		ra (float)
		dec (float)
		dir_obj (string)

		survey (str): e.g., 'hsc'
			survey of the photometric system
		z (float): 
			redshift
		pixsize (astropy angle quantity):
			in unit of arcsec

		pixelscale (astropy pixscale quantity):
			for pixel and arcsec conversion


		Decomposer Attributes
		---------------------

		bands (list): 
			e.g., ['g', 'r', 'i', 'z', 'y'] for survey = 'hsc'
		"""
		
		super(Decomposer, self).__init__(**kwargs)
		self.bands = surveysetup.surveybands[self.survey]


	def get_fp_stamp_psfmatched(self, band, bandto):
		return self.dir_obj+'stamp-{0}_psfmt-{1}.fits'.format(band, bandto)


	def get_fp_stamp_contsub(self, band, bandconti):
		return self.dir_obj+'stamp-{0}_contsub-{1}.fits'.format(band, bandconti)


	def make_stamp_linemap_I(self, bandline, bandconti, line='OIII5008', overwrite=False):
		""" 
		make stamp of line map in rest frame intensity in units of [erg s-1 cm-2 arcsec-2]
		Converted from stamp_linemap depending on self.z. Ready for isophotal measurements. 

		See make_stamp_linemap for details. 
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def make_stamp_linemap(self, bandline, bandconti, line='OIII5008', overwrite=False):
		"""
		make stamp of line map in observed frame flux in units of [erg s-1 cm-2]

		Params
		------
		self
		bandline  (str)
		bandconti (str)
		line = 'OIII5008' (str)
		overwrite = False (bool)

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., stamp-OIII5008.fits
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def make_stamp_contsub(self, band, bandconti, overwrite=True):
		"""
		make stamp that is continuum subtracted

		Params
		------
		self
		band (str)
		bandconti (str)
		overwrite=False

		Return
		------
		status

		Write Output
		------------
		e.g., stamp-i_contsub-z.fits  (if band = 'i', bandto = 'z')
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def make_stamp_psfmatch(self, band, bandto, overwrite=True):
		""" 
		make stamp that has psf matched to stamp of another band

		Params
		------
		self
		band (str)
		bandto (str)
		overwrite=False

		Return
		------
		status

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_psfmatched-z.fits
		(possibly others)
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def _get_conti_fnu_ratio_from_spector(self, band1, band2):
		""" return fnu_band1 / fnu_band2 of the continuum from spector """
		s = self._get_spector()
		ratio = s.get_fnu_ratio_band1_over_band2(band1=band1, band2=band2, component='contextrp')

		return ratio


	def _get_spector(self):
		s = spector.Spector(obj=self.obj, survey=self.survey, z=self.z)
		return s
