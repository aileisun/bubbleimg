# visual.py
# ALS 2017/06/29


import os
import numpy as np
from astropy.io import fits

import skimage.io as si


def fits_to_png(fn_in, fn_out=None, vmin=None, vmax=None, scaling='arcsinh'):
	"""
	make fits image into png image

	Params
	----------
	fn_in (string)
	fn_out = None (string)
		if set to None then it will be filenamin+'.png'
	vmin = None (float)
		the value to set to 0 (black)
	vmax = None (float)
		the value to set to 255 (white)
	scaling='arcsinh':
		linear or arcsinh


	Read in
	-------
	fn_in

	Write out
	----------
	fn_out or fn_in +'.png'
	"""

	# setting fn_out
	extension = '.png'

	if fn_out is None:
		base_in, ext_in = os.path.splitext(fn_in)

		if ext_in == '.fits': 
			fn_out = base_in+extension
		else: 
			fn_out = fn_in+extension

	if not os.path.isfile(fn_in):
		print "skipping "+fn_in+" as in file does not exist"
	else:
		# read in
		img = fits.getdata(fn_in)

		img_scaled = scale_img(img, vmin=vmin, vmax=vmax, scaling=scaling)

		si.imsave(fn_out, img_scaled)


def scale_img(img, vmin=None, vmax=None, scaling='arcsinh'):
	"""
	Params
	------
	vmin=None (float)
	vmax=None (float)
	scaling='arcsinh':
		linear or arcsinh
	"""

	epsilon = 1.e-3
	# setting vmin and vmax
	if vmin is None:
		vmin = np.min(img)
	if vmax is None:		
		vmax = np.max(img)

	# scaling
	if scaling == 'linear':
		pass
	elif scaling == 'arcsinh':
		img = np.arcsinh(img)
		vmin = np.arcsinh(vmin)
		vmax = np.arcsinh(vmax)
	else:
		raise ValueError("[visualtools] scaling not recognized")

	# clipping
	image_scaled = np.clip((img-vmin) / (vmax-vmin), 0., 1.- epsilon).astype('f2')
	image_scaled_flip = image_scaled[::-1, :]

	return image_scaled_flip
