# measurer.py 
# ALS 2017/06/29

import os
import astropy.units as u
import numpy as np
import astropy.table as at

from ..obsobj import Imager
from .. import tabtools
import noiselevel

class Measurer(Imager):

	def __init__(self, **kwargs):
		"""
		Measurer, an imager operator to do image measurement

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


		Subclass Attributes (to be overwritten by subtype)
		------------------- 
		msrtype=None (str):
			e.g., 'iso'

		"""
		
		super(Measurer, self).__init__(**kwargs)

		# subtype attribute (to be overwritten by subtype)
		self.msrtype = None


	# def get_fp_msr(self, imgtag='OIII5008_I', suffix=''):
	# 	""" return the path to the measurement results .csv file, e.g., msr_iso-OIII5008_I{suffix}.csv """
	# 	return self.dir_obj+'msr_{msrtype}-{imgtag}{suffix}.csv'.format(msrtype=self.msrtype, imgtag=imgtag, suffix=suffix)
		

	def get_fp_msr(self, msrsuffix=''):
		""" return the path to the measurement results .csv file, e.g., msr_iso.csv """
		return self.dir_obj+'msr_{msrtype}{msrsuffix}.csv'.format(msrtype=self.msrtype, msrsuffix=msrsuffix)


	# def get_fp_msrplot(self, imgtag='OIII5008_I', suffix=''):
	# 	fn_msr = self.get_fp_msr(imgtag=imgtag, suffix=suffix)
	# 	fn_noext = os.path.splitext(fn_msr)[0]
	# 	return fn_noext+'.pdf'


	def get_fp_msrplot(self, imgtag='OIII5008_I', suffix=''):
		fp_root = self.get_fp_msrtagroot(imgtag=imgtag, suffix=suffix)
		return fp_root+'.pdf'


	def get_fp_msrtagroot(self, imgtag='OIII5008_I', suffix=''):
		""" return the path to the measurement plots .pdf file, 
		e.g., msr_iso-OIII5008_I_3e-15.pdf """
		return self.dir_obj+'msr_{msrtype}-{imgtag}{suffix}'.format(msrtype=self.msrtype, imgtag=imgtag, suffix=suffix)


	# def get_fp_noiselevel(self, imgtag='OIII5008_I'):
	# 	return self.dir_obj+'noiselevel-{}.csv'.format(imgtag)


	def get_fp_noiselevel(self, msrsuffix=''):
		return self.dir_obj+'noiselevel{msrsuffix}.csv'.format(msrsuffix=msrsuffix)


	def get_fp_noiselevel_tagroot(self, imgtag='OIII5008_I'):
		return self.dir_obj+'noiselevel-{}'.format(imgtag)


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


	def make_noiselevel(self, imgtag='OIII5008_I', toplot=False, msrsuffix='', overwrite=False, append=False):
		"""
		Measure the noise level of img with tag 'imgtag' and write to noiselevel_{imgtag}.csv.

		The noise level is determined by fitting a Guassian to the histogram of the pixel values. 

		Params
		------
		self
		imgtag='OIII5008_I'
		toplot=False
		msrsuffix=''
			suffix label in the end of the measurement csv file: msr_iso_{msrsuffix}.csv.

		overwrite=False
		append=False

		Return
		------
		status
		"""
		fn = self.get_fp_noiselevel(msrsuffix=msrsuffix)
		fn_plot = self.get_fp_noiselevel_tagroot(imgtag=imgtag)+'.pdf'

		condi = {'imgtag': imgtag}

		if append or overwrite or (not tabtools.fn_has_row(fn, condi)):

			print("[measurer] making noiselevel for {}".format(imgtag))
			img = self.get_stamp_img(imgtag=imgtag, wunit=True)
			u_img = img.unit
			img = np.array(img)

			nlevel = noiselevel.getnoiselevel_gaussfit(data=img, fn_plot=fn_plot, toplot=toplot)
			tab = at.Table([[imgtag], [nlevel], [u_img.to_string()]], names=['imgtag', 'img_sigma', 'u_img'])

			tabtools.write_row(fn=fn, row=tab, condi=condi, overwrite=overwrite, append=append)
			# tab.write(fn, format='ascii.csv', overwrite=overwrite)
		else:
			print("[measurer] skip making noiselevel for {} as files exist".format(imgtag))

		return os.path.isfile(fn)


	def get_noiselevel(self, imgtag='OIII5008_I', wunit=False):
		"""
		Return the measured noise level of the image, see make_noiselevel() for details. 

		Params
		------
		self
		imgtag='OIII5008_I'
		wunit=False

		Return
		------
		n_level (float or quantity)
		"""
		fn = self.get_fp_noiselevel()

		self.make_noiselevel(imgtag=imgtag, toplot=False, overwrite=False)

		tab = at.Table.read(fn, format='ascii.csv')
		tab[tab['imgtag'] == imgtag]

		n_level = tab['img_sigma'][0]

		if wunit:
			u_img = tab['u_img'][0]
			n_level = n_level * u.Unit(u_img)

		return n_level


