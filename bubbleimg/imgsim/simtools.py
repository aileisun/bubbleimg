# psf.py

import numpy as np
from skimage.util import random_noise
from astropy.io import fits
# from bubblepy.measureimg import gaussfit
# from scipy.ndimage.filters import convolve


def add_gaussnoise(img, noise_sigma):
	""" 
	add iid gaussian noise on image

	Params
	------
	img
	noise_sigma

	Return
	------
	img
	"""
	if noise_sigma == 0:
		result = img
	else:
		result = img+random_noise(np.zeros(img.shape), mode='gaussian',var=noise_sigma**2,clip=False)

	return result

def bin_img(img, binsize=2):
	""" 
	bin the image by taking the mean average of (binsize * binsize) blocks. The resulting image is binsize times smaller than the original in each dimension. If the img size is not divisible by the binsize, the remainder is cut off and ignored. 

	Params
	------
	img
	binsize = 2 (int)

	Return
	------
	img
	"""

	nx, ny = img.shape

	# sanity check
	for n in [nx, ny]:
		if n%binsize != 0:
			print("[simtools] image size not divisible by binsize. The remainder is ignored. ")
			# 		raise Exception("[simtools] image size is not divisible by binsize")

	nx_b, ny_b = nx // binsize, ny // binsize

	img_b = np.zeros([nx_b, ny_b])

	for i in range(nx_b):
		for j in range(ny_b):
			img_b[i, j] = np.average(img[i*binsize:(i+1)*binsize, j*binsize:(j+1)*binsize])

	return img_b

