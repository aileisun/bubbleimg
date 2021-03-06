# imager.py 
# ALS 2017/10/05

import os
import astropy.units as u
from astropy.io import fits
import numpy as np

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from .operator import Operator
from .. import standards
from ..filters import surveysetup
from ..external_links import file_humvi_compose
from .. import visualtools


class Imager(Operator):

	def __init__(self, **kwargs):
		"""
		Imager, parent class for all obj operator for images

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

		z=-1 (float):  
			redshift, if not provided, use self.obj.z or self.obj.sdss.z. It does not automatically query sdss to get z. If nothing is pecified then set to -1. 

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

		"""
		
		super(Imager, self).__init__(**kwargs)

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
		elif 'z' in kwargs:
			self.z = kwargs.pop('z') 
		else: 
			print("[imager] not redshift used, assuming -1")
			self.z = -1

		# set center_mode
		self.center_mode = kwargs.pop('center_mode', 'n/2')

		# set pixsize
		self.pixsize = surveysetup.pixsize[self.survey]
		self.pixelscale = u.pixel_scale(self.pixsize/u.pixel)


	def get_fp_stamp(self, band):
		return self.dir_obj + self.get_fn_stamp(band)


	def get_fn_stamp(self, band):
		return 'stamp-{0}.fits'.format(band)


	def get_fp_psf(self, band):
		return self.dir_obj + self.get_fn_psf(band)


	def get_fn_psf(self, band):
		return 'psf-{0}.fits'.format(band)


	def get_fp_stamp_line(self, line):
		""" e.g., stamp-OIII5008.fits, for stamp in observed frame in flux """
		return self.dir_obj+'stamp-{0}.fits'.format(line)


	def get_fp_stamp_line_I(self, line):
		""" e.g., stamp-OIII5008_I.fits for stamp in rest frame in intensity"""
		return self.dir_obj+'stamp-{0}_I.fits'.format(line)


	def get_fp_stamp_img(self, imgtag):
		""" e.g., stamp-{imgtag}.fits"""
		return self.dir_obj+'stamp-{imgtag}.fits'.format(imgtag=imgtag)


	def get_stamp_img(self, imgtag, wunit=False):
		""" return the image as numpy array """
		fn_img = self.get_fp_stamp_img(imgtag=imgtag)
		hdus = fits.open(fn_img)
		if wunit:
			img = hdus[0].data * u.Unit(hdus[0].header['BUNIT'])
		else: 
			img = hdus[0].data

		return img


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


	def make_colorimg(self, bands ='riz', img_type='stamp', s=[1.0, 1.1, 1.0], p=[1.6, 1.6], overwrite=False):
		"""
		make color composit image using external package HumVI. Example file name: 'color_stamp-riz.png'.

		See Humvi documentation -- https://github.com/drphilmarshall/HumVI

		Params
		------
		bands ='riz'
		img_type='stamp'
		s=[1.0, 1.1, 1.0]
			humvi params for color balance. Default is set to bring out the middle band. 
		p=[1.6, 1.6]
			humvi params Q, alpha
		overwrite=False

		Return
		------
		status (bool)
		"""
		fn = self.dir_obj+'color_{img_type}-{bands}.png'.format(bands=bands, img_type=img_type)

		fns_in = [self.dir_obj+img_type+'-'+band+'.fits' for band in bands[::-1]]

		# set params, for example, commandparams = ' -s 1.0,1.1,1.0  -p 1.6,1.6  -o '
		commandparams = ' -s {}  -p {}  -o '.format(','.join(np.array(s).astype('str')), ','.join(np.array(p).astype('str')))

		if (not os.path.isfile(fn)) or overwrite:
			commandfiles = '{0} {1} {2} {3}'.format(fn, fns_in[0], fns_in[1], fns_in[2])
			commandHumVI = file_humvi_compose+commandparams+commandfiles

			os.system(commandHumVI)

		status = os.path.isfile(fn)

		return status


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
