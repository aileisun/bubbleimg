# simulator.py 
# ALS 2017/10/05
import os
from astropy.io import fits

from ..obsobj import Imager
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


	def get_fp_noised(self, imgtag, img_sigma):
		return self.get_fp_stamp_img(imgtag+'_noised-{}'.format('%.1f'%img_sigma))


	def make_noised(self, imgtag, img_sigma=1, overwrite=False):
		fn = self.get_fp_noised(imgtag=imgtag, img_sigma=img_sigma)

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

