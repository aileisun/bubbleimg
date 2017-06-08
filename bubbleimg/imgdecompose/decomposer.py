# decomposer.py 
# ALS 2017/06/01

from ..obsobj import Operator
from .. import filters
from .. import spector
from .. import imgdownload


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
		self.bands = filters.filtertools.surveybands[self.survey]


	def get_fp_stamp(self, band):
		L = imgdownload.imgLoader(obj=self.obj)
		return L.get_fp_stamp(band)
		# return self.dir_obj+'stamp-{0}.fits'.format(band)


	def get_fp_psf(self, band):
		L = imgdownload.imgLoader(obj=self.obj)
		return L.get_fp_psf(band)
		# return self.dir_obj+'psf-{0}.fits'.format(band)


	def get_fp_stamp_psfmatched(self, band, bandto):
		return self.dir_obj+'stamp-{0}_psfmatched-{1}.fits'.format(band, bandto)


	def get_fp_stamp_contsub(self, band, bandconti):
		return self.dir_obj+'stamp-{0}_contsub-{1}.fits'.format(band, bandconti)


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

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_contsub-z.fits
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def _get_conti_fnu_ratio_from_spector(self, band1, band2):
		""" return fnu_band1 / fnu_band2 of the continuum from spector """
		s = spector.Spector(obj=self.obj, survey=self.survey, z=self.z)
		ratio = s.get_fnu_ratio_band1_over_band2(band1=band1, band2=band2, component='contextrp')

		return ratio


