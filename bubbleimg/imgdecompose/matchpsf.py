# matchpsf.py
# ALS 2017/06/05

import os
import astropy.convolution as ac
import scipy.optimize as so
import numpy as np

def match_psf(img, psf, psfto): 
	"""
	match the psf of the image from psf to psfto

	Params
	------
	img
	psf
	psfto

	Return
	------
	img_cnvled:
		the psf matched image
	psf_cnvled:
		the matched psf (psf convolved with kernel_cnvl)
	kernel_cnvl:
		the moffat kernel that can smear psf into psfto
	"""

	# get convolving kernel
	kernel_cnvl = get_convolving_moffat_kernel(psf, psfto)

	# make psf matched img
	img_cnvled = ac.convolve(img, kernel_cnvl) 

	# make psf matched psf
	psf_cnvled = ac.convolve(psf, kernel_cnvl)

	return img_cnvled, psf_cnvled, kernel_cnvl



def get_convolving_moffat_kernel(arr_in, arr_out, nx=43, ny=43, scale_factor=10000., tolerance=1.e-3):
	"""
	find k such that convolve(arr_in, k) = arr_out

	Params
	------
	arr_in
	arr_out

	Return
	------
	k:
		convolving kernel
	"""

	def cost_func(plist):
		"""
		chisq of kernel_to - convl(kernel_from, kernel_via) in kernel space.
		"""
		gamma, alpha = plist
		k = ac.Moffat2DKernel(gamma, alpha, x_size=nx, y_size=ny)

		arr_out_predict = ac.convolve(arr_in, k)

		arr_out_fit, arr_out_predict_fit = match_dimension(arr_out, arr_out_predict)
		diff = (arr_out_fit - arr_out_predict_fit)*scale_factor

		return np.sum(diff**2)/diff.size

	bounds = ((1.e-2, np.inf), (-np.inf, np.inf))
	plist_init = (3., 1.)
	res = so.minimize(cost_func, plist_init, bounds=bounds)


	gamma_via, alpha_via = res.x

	avg_diff = np.sqrt(res.fun)/scale_factor
	if avg_diff > tolerance:
		warnings.warn("[psfmatch.deconvl_moffat_fast] best fit worse than tolerance. ", UserWarning)

	k = ac.Moffat2DKernel(gamma_via, alpha_via, x_size=nx, y_size=ny)

	return k


def isodd(x):
	return x % 2 != 0


def match_dimension(p0, p1):
	"""
	Params
	------
	p0 (array)
	p1 (array)

	Return
	------
	p0 (array)
	p1 (array)
	"""

	if p0.shape != p1.shape:
		nxmax = max(p0.shape[0], p1.shape[0])
		nymax = max(p0.shape[1], p1.shape[1])

		p0 = pad_edge_to_shape(p0, nxmax, nymax)
		p1 = pad_edge_to_shape(p1, nxmax, nymax)

	return p0, p1


def pad_edge_to_shape(arr, nx, ny):
	""" 
	pad edge of array such that its dimension becomes nx, ny. 
	only works for odd number length arrays
	"""
	nxa, nya = arr.shape

	for n in [nx, ny, nxa, nya]:
		if not isodd(n):
			raise Exception("[pad_edge_to_shape] only accept odd-number-length arrays")

	if nx > nxa:
		dx = (nx - nxa)/2
		arr = np.pad(arr, ((dx, dx), (0, 0)), 'constant', constant_values=0)
	elif nx < nxa:
		raise Exception("[pad_edge_to_shape] final dimension smaller than array")

	if ny > nya:
		dy = (ny - nya)/2
		arr = np.pad(arr, ((0, 0), (dy, dy)), 'constant', constant_values=0)
	elif ny < nya:
		raise Exception("[pad_edge_to_shape] final dimension smaller than array")

	return arr


def normalize_kernel(kernel):
	return kernel/np.sum(kernel)


