# simulator.py 
# ALS 2017/10/05
import os
from astropy.io import fits
import numpy as np

from ..obsobj import Imager
from .. import imgmeasure
from .. import tabtools
import psfnoisetools

class Simulator(Imager):

	def __init__(self, **kwargs):
		"""
		Simulator, an imager operator to do image simulation, such as adding noise

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

		"""
		
		super(Simulator, self).__init__(**kwargs)


	def get_tag_noised(self, img_sigma):
		""" e.g., '_noised-1.0' """
		return '_noised-{}'.format('%.1f'%img_sigma)


 	def get_fp_noised(self, imgtag='OIII5008_I', img_sigma=1, suffix=''):
 		""" e.g., 'stamp-OIII5008_I_noised-1.0.fits' """
		return self.get_fp_stamp_img(imgtag+self.get_tag_noised(img_sigma)+suffix)


	def _get_measurer(self, msrtype='iso'):
		if msrtype=='iso':
			return imgmeasure.isoMeasurer(obj=self.obj, survey=self.survey, z=self.z, center_mode=self.center_mode)
		else: 
			raise InputError("[simulator] msrtype not understood")


	def make_noised(self, imgtag='OIII5008_I', img_sigma=1, suffix='', overwrite=False):
		"""
		create noised image, e.g., stamp-OIII5008_I_noised-1.0.fits

		Params
		------
		imgtag='OIII5008_I' (str)
		img_sigma=1 (float)
		suffix='' (str)
			to be attached to the end of the output filename
		overwrite=False

		Return
		------
		status (bool)

		Write Output
		------------
		e.g., stamp-OIII5008_I_noised-1.0.fits
		"""
		fn = self.get_fp_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=suffix)

		if not os.path.isfile(fn) or overwrite:
			fn_img = self.get_fp_stamp_img(imgtag=imgtag)
			hdus = fits.open(fn_img)
			img = hdus[0].data
			img = psfnoisetools.add_gaussnoise(img=img, noise_sigma=img_sigma)

			# write
			hdus[0].data = img
			hdus[0].header['COMMENT'] = 'added noise of sigma = {}'.format('%.2f'%img_sigma)
			hdus.writeto(fn, overwrite=overwrite)

		else:
			print("[simulator] skip making noise")

		return os.path.isfile(fn)


	def sim_noised(self, imgtag='OIII5008_I', img_sigma=1, niter=100, msrtype='iso', running_indx=False, keep_img=False, summarize=True, overwrite=False, **msrkwargs):
		"""
		simulate noise images and make measurements with 'niter' iterations. 

		Params
		------
		imgtag='OIII5008_I' (str)
		img_sigma=1 (float)
		niter=100
		msrtype='iso'
		running_indx=False: 
			give each of the noised image .fits file an index, starting from 0. 
		keep_img=False
			whether to keep fits images
		summarize=True
			produce summarizing file as _smr.csv files. 
		overwrite=False
		**msrkwargs:
			additional arguments for m.make_measurements()

		Return
		------
		status (bool)
 
		Write Output
		------------
		e.g., msr_iso_noised-1.0.csv, noiselevel_noised-1.0.csv
			\if summarize, then also
			msr_iso_noised-1.0_smr.csv, noiselevel_noised-1.0_smr.csv
		"""
		m = self._get_measurer(msrtype=msrtype)
		tag_noised = self.get_tag_noised(img_sigma=img_sigma)
		fn_msr = m.get_fp_msr(msrsuffix=tag_noised)
		fn_nl = m.get_fp_noiselevel(msrsuffix=tag_noised)
		img_suffix = ''

		if not os.path.isfile(fn_msr) or overwrite:
			if (os.path.isfile(fn_msr)) and overwrite:
				# Delete before writing
				os.remove(fn_msr)
				os.remove(fn_nl)

			for i in range(niter):
				if running_indx:
					img_suffix = '_{}'.format(str(i))

				fn_img = self.get_fp_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=img_suffix)
				status = self.make_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=img_suffix, overwrite=True)

				if status:
					m.make_noiselevel(imgtag=imgtag+tag_noised+img_suffix, toplot=False, msrsuffix=tag_noised, overwrite=False, append=True)

					m.make_measurements(imgtag=imgtag+tag_noised+img_suffix, savecontours=False, plotmsr=False, msrsuffix=tag_noised, overwrite=False, append=True, **msrkwargs)

					if summarize:
						status = self._summarize_sim_noised_kernel(imgtag, tag_noised, img_suffix, **msrkwargs)
				else:
					raise Exception("[simulator] noised image not successfully created")

				if not keep_img:
					os.remove(fn_img)

		else: 
			print("[simulator] skip sim_noised as file exists")

		return os.path.isfile(fn_msr)


	def _summarize_sim_noised_kernel(self, imgtag, tag_noised, img_suffix, **msrkwargs):
		m = self._get_measurer()
		imgtagkwargs = {'imgtag': imgtag+tag_noised+img_suffix}
		fn_msr = m.get_fp_msr(msrsuffix=tag_noised)
		fn_nl = m.get_fp_noiselevel(msrsuffix=tag_noised)
		fn_msr_smr = m.get_fp_msr_smr(msrsuffix=tag_noised)
		fn_nl_smr = m.get_fp_noiselevel_smr(msrsuffix=tag_noised)
		img_suffix = ''

		condi = dict(msrkwargs.items() + imgtagkwargs.items())
		status1 = tabtools.summarize(fn_in=fn_msr, fn_out=fn_msr_smr, columns=[], condi=condi, overwrite=True)

		condi = imgtagkwargs
		status2 = tabtools.summarize(fn_in=fn_nl, fn_out=fn_nl_smr, columns=[], condi=condi, overwrite=True)

		return np.all([status1, status2])


	def summarize_sim_noised(self, imgtag='OIII5008_I', img_sigma=1, msrtype='iso', overwrite=False, **msrkwargs):
		"""
		stand alone function to summarize measurements and noiselevels. 

		Params
		------
		img_sigma=1 (float)
		msrtype='iso'
		columns=[]
			list of column names to summarize over
		overwrite=False
		"""
		m = self._get_measurer(msrtype=msrtype)
		tag_noised = self.get_tag_noised(img_sigma=img_sigma)
		img_suffix = ''

		status = self._summarize_sim_noised_kernel(imgtag, tag_noised, img_suffix, **msrkwargs)
		return status

