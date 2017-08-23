# measurer.py 
# ALS 2017/06/29

import os
import astropy.units as u

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from ..obsobj import Operator
from .. import imgdecompose
from .. import standards
from ..filters import surveysetup
from ..imgdownload.loader import imgLoader


class Measurer(Operator):

	def __init__(self, **kwargs):
		"""
		Measurer, an obj operator to do image measurement

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

		center_mode='n/2' (str):
			how is image center defined in case of even n, 'n/2' or 'n/2-1'. Should be set to n/2-1 if the image is downloaded from HSC quarry. 


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

		pixelscale (astropy pixscale quantity):
			for pixel and arcsec conversion

		Subclass Attributes: (to be overwritten by subtype)
			msrtype=None (str):
				e.g., 'iso'

		"""
		
		super(Measurer, self).__init__(**kwargs)

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

		# set center_mode
		self.center_mode = kwargs.pop('center_mode', 'n/2')

		# set pixsize
		self.pixsize = surveysetup.pixsize[self.survey]
		self.pixelscale = u.pixel_scale(self.pixsize/u.pixel)

		# subtype attribute (to be overwritten by subtype)
		self.msrtype = None


	def get_fp_stamp_line(self, line):
		""" e.g., stamp-OIII5008.fits, for stamp in observed frame in flux """
		d = self._get_decomposer()
		return d.get_fp_stamp_line(line)


	def get_fp_stamp_line_I(self, line):
		""" e.g., stamp-OIII5008_I.fits for stamp in rest frame in intensity"""
		d = self._get_decomposer()
		return d.get_fp_stamp_line_I(line)


	def get_fp_stamp_img(self, imgtag):
		""" e.g., stamp-{imgtag}.fits"""
		return self.dir_obj+'stamp-{imgtag}.fits'.format(imgtag=imgtag)


	def get_fp_msr(self, imgtag='OIII5008_I', suffix=''):
		""" return the path to the measurement results .csv file, e.g., msr_iso-OIII5008_I{suffix}.csv """
		return self.dir_obj+'msr_{msrtype}-{imgtag}{suffix}.csv'.format(msrtype=self.msrtype, imgtag=imgtag, suffix=suffix)
		

	def get_fp_msrplot(self, imgtag='OIII5008_I', suffix=''):
		fn_msr = self.get_fp_msr(imgtag=imgtag, suffix=suffix)
		fn_noext = os.path.splitext(fn_msr)[0]
		return fn_noext+'.pdf'


	def make_measurements_line_I(self, line='OIII5008', overwrite=False, **kwargs):
		"""
		make measurements on the intensity map of a line. 
			e.g., stamp-OIII5008_I.fits
		See make_measurements for details
		"""
		self.make_measurements(imgtag=line+'_I', overwrite=overwrite, **kwargs)


	def make_measurements(self, imgtag='OIII5008_I', overwrite=False, **kwargs):
		"""
		make measurements on a map
			if imgtag='OIII5008_I' then measure 'stamp-OIII5008_I.fits'

		Params
		------
		self
		imgtag='OIII5008_I'
		overwrite = False (bool)

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., msr_iso-OIII5008.csv
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	def make_colorimg(self, bands ='riz', img_type='stamp', overwrite=False):
		"""
		make color composit image using external package HumVI. Example file name: 'color_stamp-riz.png'.
		Uses class method of loader. 

		Params
		------
		bands ='riz'
		img_type='stamp'
		overwrite=False

		Return
		------
		status (bool)
		"""
		l = imgLoader(obj=self.obj)
		status = l.plot_colorimg(bands=bands, img_type=img_type, overwrite=overwrite)

		return status


	def _get_decomposer(self):
		d = imgdecompose.Decomposer(obj=self.obj, survey=self.survey, z=self.z)
		return d


	def _theta_to_pix(self, theta):
		""" convert an angular quantity (*u.arcsec) theta to pixel scale (float) """
		return (theta.to(u.pix, self.pixelscale)/u.pix).to(u.dimensionless_unscaled)


	def _pix_to_theta(self, pix, wunit=True):
		""" convert pix (float) to angular quantity (*u.arcsec) """
		result = (pix*u.pix).to(u.arcsec, self.pixelscale)
		if wunit:
			return result
		else: 	
			return (result/u.arcsec).to(u.dimensionless_unscaled)


	def _theta_to_kpc(self, theta, wunit=True):
		result = (theta*self._get_kpc_proper_per_arcsec()).to(u.kpc)

		if wunit:
			return result
		else: 	
			return (result/u.kpc).to(u.dimensionless_unscaled)


	def _get_kpc_proper_per_arcsec(self):
		return 1./cosmo.arcsec_per_kpc_proper(self.z)


	def _get_xc_yc(self, img):
		""" return (xc, yc) the coordinate of image center given image, according to self.center_mode """

		xc, yc = standards.get_img_xycenter(img, center_mode=self.center_mode)
		return xc, yc