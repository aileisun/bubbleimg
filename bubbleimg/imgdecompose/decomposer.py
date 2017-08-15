# decomposer.py 
# ALS 2017/06/01

import os

from ..obsobj import Operator
from .. import spector
from .. import imgdownload
from ..filters import surveysetup
from .. import visualtools


class Decomposer(Operator):

	def __init__(self, **kwargs):
		"""
		Decomposer, an obj operator to do image decomposition

		Params
		------
		Operator params:
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
			redshift, if not provided, use self.obj.z or self.obj.sdss.z. It does not automatically query sdss to get z. 


		Attributes
		----------
		Operator Attributes:	
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
		bands (list): 
			e.g., ['g', 'r', 'i', 'z', 'y'] for survey = 'hsc'
		"""
		
		super(Decomposer, self).__init__(**kwargs)

		# set survey
		if hasattr(self.obj, 'survey'):
			default_survey = self.obj.survey
			self.survey = kwargs.pop('survey', default_survey)
		else: 
			self.survey = kwargs.pop('survey')

		# set up obj.survey
		if self.survey == 'hsc':
			self.obj.add_hsc()
		elif self.survey == 'sdss':
			self.obj.add_sdss()

		# set z
		if hasattr(self.obj, 'z'):
			self.z = kwargs.pop('z', self.obj.z)
		elif hasattr(self.obj, 'sdss'):
			self.z = kwargs.pop('z', self.obj.sdss.z) 
		else: 
			self.z = kwargs.pop('z') 

		# set other attributes
		self.pixsize = surveysetup.pixsize[self.survey]
		self.bands = surveysetup.surveybands[self.survey]



	def get_fp_stamp(self, band):
		L = imgdownload.imgLoader(obj=self.obj)
		return L.get_fp_stamp(band)


	def get_fp_stamp_psfmatched(self, band, bandto):
		return self.dir_obj+'stamp-{0}_psfmt-{1}.fits'.format(band, bandto)


	def get_fp_stamp_contsub(self, band, bandconti):
		return self.dir_obj+'stamp-{0}_contsub-{1}.fits'.format(band, bandconti)


	def get_fp_stamp_line(self, line):
		""" e.g., stamp-OIII5008.fits, for stamp in observed frame in flux """
		return self.dir_obj+'stamp-{0}.fits'.format(line)


	def get_fp_stamp_line_I(self, line):
		""" e.g., stamp-OIII5008_I.fits for stamp in rest frame in intensity"""
		return self.dir_obj+'stamp-{0}_I.fits'.format(line)


	def plot_stamp_linemap_I(self, line='OIII5008', overwrite=False, vmin=None, vmax=10.):
		""" 
		plot line map I as png. 

		Params
		------
		self
		line='OIII5008' (str)
		overwrite=False (bool)

		Return
		------
		status (bool)
		"""
		fn = self.get_fp_stamp_line_I(line=line)
		fn_out = os.path.splitext(fn)[0]+'.png'

		if not os.path.isfile(fn_out) or overwrite:
			print("[decomposer] plotting linemap")
			visualtools.fits_to_png(fn_in=fn, fn_out=fn_out, vmin=vmin, vmax=vmax, scaling='arcsinh')
		else: 
			print("[decomposer] skip plotting linemap as files exist")

		return os.path.isfile(fn_out)


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
