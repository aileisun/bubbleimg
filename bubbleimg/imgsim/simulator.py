# simulator.py 
# ALS 2017/10/05
import os
from astropy.io import fits
import numpy as np
import astropy.convolution as ac
import astropy.table as at
import copy

from ..obsobj import Imager
from .. import imgmeasure
from .. import imgdecompose
from ..imgdecompose.plain import matchpsf
from .. import tabtools
import simtools

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


	def get_fp_stamp_noised(self, imgtag='OIII5008_I', img_sigma=1, suffix=''):
 		""" e.g., 'stamp-OIII5008_I_noised-1.0.fits' """
		return self.get_fp_stamp_img(imgtag+self.get_tag_noised(img_sigma)+suffix)


	def get_fp_stamp_binned(self, imgtag='OIII5008_I', binsize=2):
		return self.get_fp_stamp_img(imgtag+'_binned-{binsize}'.format(binsize=str(binsize)))


	def get_fp_psf_binned(self, imgtag='OIII5008_I', binsize=2):
		return self.get_fp_psf(imgtag+'_binned-{binsize}'.format(binsize=str(binsize)))


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
		fn = self.get_fp_stamp_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=suffix)

		if not os.path.isfile(fn) or overwrite:
			fn_img = self.get_fp_stamp_img(imgtag=imgtag)
			hdus = fits.open(fn_img)
			img = hdus[0].data
			img = simtools.add_gaussnoise(img=img, noise_sigma=img_sigma)

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

				fn_img = self.get_fp_stamp_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=img_suffix)
				status = self.make_noised(imgtag=imgtag, img_sigma=img_sigma, suffix=img_suffix, overwrite=True)

				if status:
					m.make_noiselevel(imgtag=imgtag+tag_noised+img_suffix, toplot=False, msrsuffix=tag_noised, overwrite=False, append=True)

					m.make_measurements(imgtag=imgtag+tag_noised+img_suffix, savecontours=False, plotmsr=False, msrsuffix=tag_noised, overwrite=False, append=True, **msrkwargs)

					if summarize:
						status = self._summarize_sim_noised_kernel(imgtag, tag_noised, img_suffix, overwrite=True, **msrkwargs)
				else:
					raise Exception("[simulator] noised image not successfully created")

				if not keep_img:
					os.remove(fn_img)

		else: 
			print("[simulator] skip sim_noised as file exists")

		return os.path.isfile(fn_msr)


	def _summarize_sim_noised_kernel(self, imgtag, tag_noised, img_suffix, overwrite=False, **msrkwargs):
		m = self._get_measurer()
		imgtagkwargs = {'imgtag': imgtag+tag_noised+img_suffix}
		fn_msr = m.get_fp_msr(msrsuffix=tag_noised)
		fn_nl = m.get_fp_noiselevel(msrsuffix=tag_noised)
		fn_msr_smr = m.get_fp_msr_smr(msrsuffix=tag_noised)
		fn_nl_smr = m.get_fp_noiselevel_smr(msrsuffix=tag_noised)
		img_suffix = ''

		condi = dict(msrkwargs.items() + imgtagkwargs.items())
		status1 = tabtools.summarize(fn_in=fn_msr, fn_out=fn_msr_smr, columns=[], condi=condi, overwrite=overwrite)

		condi = imgtagkwargs
		status2 = tabtools.summarize(fn_in=fn_nl, fn_out=fn_nl_smr, columns=[], condi=condi, overwrite=overwrite)

		return np.all([status1, status2])


	def summarize_sim_noised(self, imgtag='OIII5008_I', img_sigma=1, msrtype='iso', overwrite=False, **msrkwargs):
		"""
		stand alone function to summarize measurements and noiselevels. Always overwrites. 

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

		status = self._summarize_sim_noised_kernel(imgtag, tag_noised, img_suffix, overwrite=overwrite, **msrkwargs)
		return status


	def make_binned(self, imgtag='OIII5008_I', binsize=2, binpsf=False, overwrite=False):
		"""
		create binned image, e.g., stamp-OIII5008_I_binned-2.fits

		Params
		------
		imgtag='OIII5008_I' (str)
		binsize = 2
			how many pixels to bin into one pix
		binpsf=False
			whether to also produce binned psf
		overwrite=False

		Return
		------
		status (bool)

		Write Output
		------------
		e.g., stamp-OIII5008_I_binned-2.fits
		"""
		fn = self.get_fp_stamp_binned(imgtag=imgtag, binsize=binsize)

		if not os.path.isfile(fn) or overwrite:
			fn_img_in = self.get_fp_stamp_img(imgtag=imgtag)
			
			hdus = fits.open(fn_img_in)
			img = hdus[0].data
			img = simtools.bin_img(img=img, binsize=binsize)
			# write
			hdus[0].data = img
			hdus[0].header['COMMENT'] = 'binned with binsize = {}'.format(str(binsize))
			hdus.writeto(fn, overwrite=overwrite)

			#=== psf
			if binpsf:
				print("[simulator] WARNING: bin psf not tested")
				fn_psf = self.get_fp_psf_binned(imgtag=imgtag, binsize=binsize)
				fn_psf_in = self.get_fp_psf(band=imgtag)
				hdus = fits.open(fn_psf_in)
				img = hdus[0].data
				img = simtools.bin_img(img=img, binsize=binsize)
				# write
				hdus[0].data = img
				hdus[0].header['COMMENT'] = 'binned with binsize = {}'.format(str(binsize))
				hdus.writeto(fn_psf, overwrite=overwrite)

		else:
			print("[simulator] skip making binned image a file exists")

		return os.path.isfile(fn)


	def get_tag_smeared(self, gamma=3., alpha=1.):
		return '_smeared-moffatg{gamma}a{alpha}'.format(gamma='%.0f'%gamma, alpha='%.0f'%alpha,)


	def make_smeared(self, imgtag='OIII5008_I', gamma=3., alpha=1., nx=43, ny=43, overwrite=True):
		"""
		make image convolved by moffat psf kernel. 

		Params
		------
		self, imgtag='OIII5008_I', gamma=3., alpha=1., nx=43, ny=43, overwrite=True

		Return
		------
		status (bool)

		Write Output
		------------
		e.g., 
			stamp-OIII5008_I_smeared-moffatg2a1.fits
			psf-OIII5008_I_smeared-moffatg2a1.fits
		"""
		tag_smeared = self.get_tag_smeared(gamma=gamma, alpha=alpha)
		fn = self.get_fp_stamp_img(imgtag=imgtag+tag_smeared)

		if not os.path.isfile(fn) or overwrite:
			d = self._get_decomposer()

			fn_psf = d.get_fp_psf(fn)
			fn_stamp_in = self.get_fp_stamp_img(imgtag=imgtag)
			fn_psf_in = d.get_fp_psf(fn_stamp_in)

			img = fits.getdata(fn_stamp_in)
			psf = fits.getdata(fn_psf_in)

			kernel = ac.Moffat2DKernel(gamma, alpha, x_size=nx, y_size=ny)
			img_smeared = ac.convolve(img, kernel)
			psf_smeared = ac.convolve(psf, kernel)

			comment = 'PSF smeared by moffat kernel with gamma {} alpha {}'.format(str(gamma), str(alpha))
			matchpsf.replace_img_in_fits(fn_stamp_in, fn, img_smeared, comment=comment, overwrite=overwrite)
			fits.PrimaryHDU(psf_smeared).writeto(fn_psf, overwrite=overwrite)

		else:
			print("[simulator] skip making smeared image a file exists")

		return os.path.isfile(fn)


	def sim_smeared(self, imgtag='OIII5008_I', smearargs=at.Table([[2.], [2.]], names=['gamma', 'alpha']), msrtype='iso', keep_img=True, overwrite=True, **msrkwargs):
		"""
		simulate noise images and make measurements with 'niter' iterations. 

		Params
		------
		imgtag='OIII5008_I' (str)
		smearargs=at.Table([[1.], [1.]], names=['gamma', 'alpha']) (astropy table)
			a table containing the param sets to use for smearing
		msrtype='iso'
		keep_img=False
			whether to keep fits images
		overwrite=False
		**msrkwargs:
			additional arguments for m.make_measurements()

		Return
		------
		status (bool)
 
		Write Output
		------------
		e.g., msr_iso_smeared-moffatg1a1.csv, psf_smeared-moffatg1a1.csv
		"""

		m = self._get_measurer(msrtype=msrtype)
		d = self._get_decomposer()
		msrsuffix = '_smeared'
		fn_msr = m.get_fp_msr(msrsuffix=msrsuffix)
		fn_psftab = d.get_fp_psf_tab(msrsuffix=msrsuffix)

		if not os.path.isfile(fn_msr) or overwrite:
			if (os.path.isfile(fn_msr)) and overwrite:
				# Delete before writing
				os.remove(fn_msr)
				os.remove(fn_psftab)

			# measure before smearing
			m.make_measurements(imgtag=imgtag, savecontours=False, plotmsr=False, msrsuffix=msrsuffix, overwrite=False, append=True, **msrkwargs)
			d.make_psf_tab(imgtag=imgtag, msrsuffix=msrsuffix, overwrite=False, append=True)

			for row in smearargs:
				gamma = row['gamma']
				alpha = row['alpha']

				tag_smeared = self.get_tag_smeared(gamma=gamma, alpha=alpha)
				status = self.make_smeared(imgtag=imgtag, gamma=gamma, alpha=alpha, overwrite=True)

				if status:
					m.make_measurements(imgtag=imgtag+tag_smeared, savecontours=False, plotmsr=False, msrsuffix=msrsuffix, overwrite=False, append=True, **msrkwargs)
					d.make_psf_tab(imgtag=imgtag+tag_smeared, msrsuffix=msrsuffix, overwrite=False, append=True)
				else:
					raise Exception("[simulator] smeared image not successfully created")

				if not keep_img:
					fn_img = self.get_fp_stamp_img(imgtag=imgtag+tag_smeared)
					fn_psf = d.get_fp_psf(fn_img)
					os.remove(fn_img)
					os.remove(fn_psf)

			# add header
			theader_0 = at.Table([[0.], [0.]], names=['gamma', 'alpha'])
			theader = at.vstack([theader_0, smearargs])
			theader.rename_column('gamma', 'smearing_gamma')
			theader.rename_column('alpha', 'smearing_alpha')
			_attach_header_columns_to_tab(fn_msr, theader)
			_attach_header_columns_to_tab(fn_psftab, theader)

		else: 
			print("[simulator] skip sim_smeared as file exists")

		return os.path.isfile(fn_msr)


	def _get_measurer(self, msrtype='iso'):
		if msrtype == 'iso':
			return imgmeasure.isoMeasurer(obj=self.obj, survey=self.survey, z=self.z, center_mode=self.center_mode)
		else: 
			raise InputError("[simulator] msrtype not understood")


	def _get_decomposer(self, decomtype='plain'):
		if decomtype == 'plain':
			return imgdecompose.plainDecomposer(obj=self.obj, survey=self.survey, z=self.z, center_mode=self.center_mode)
		else: 
			raise InputError("[simulator] decomtype not understood")

def _attach_header_columns_to_tab(fn, theader):
	tab = at.Table.read(fn)
	tab = at.hstack([theader, tab])
	tab.write(fn, overwrite=True)
