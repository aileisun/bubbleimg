# matchpsf.py
# ALS 2017/06/05

import os
import numpy as np

import scipy.optimize as so
import scipy.interpolate as si

import astropy.convolution as ac
import astropy.modeling as am
from astropy.io import fits


def match_psf_fits(fp_img, fp_psf, fp_psfto, fp_img_out, fp_psf_out, fp_psk_out, overwrite=True, towrite_psk=False):
	""" 
	match the psf of the image to that of psfto, given fits file names. For details, see match_psf().

	Params
	------
	fp_img		(str):
		(read) filepath to the image to be psf matched
	fp_psf		(str):
		(read) filepath to the psf of image
	fp_psfto	(str):
		(read) filepath to the psf to be matched to
	fp_img_out	(str):
		(write) filepath to the psf matched image output
	fp_psf_out	(str):
		(write) filepath to the psf of the matched image
	fp_psk_out	(str):
		(write) filepath to the psf matching kernel
	overwrite=True (bool)
	towrite_psf=False (bool):
		whether to write psf matching kernel

	Write output
	------------
	fp_img_out
	fp_psf_out
	fp_psk_out
	"""

	img = fits.getdata(fp_img)
	psf = fits.getdata(fp_psf)
	psfto = fits.getdata(fp_psfto)

	# sanity check -- normalization of psf
	for a in [psf, psfto]:
		if np.absolute(np.sum(a) - 1.) > 3.e-4:
			raise ValueError("[matchpsf] input psf is not normalized with sum".format('%.5f'%np.sum(a)))

	img_out, psf_out, psk_out, = match_psf(img, psf, psfto)

	replace_img_in_fits(fn_from=fp_img, fn_to=fp_img_out, img=img_out, comment="PSF matched by ALS", overwrite=overwrite)

	# # construct hdus for outputs
	# hdus = fits.open(fp_img)
	# hdus[0].data = img_out
	# hdus[0].header['COMMENT'] = "PSF matched by ALS"
	# # hdus[0].header['COMMENT'] = "    from {} to {}".format(os.path.basename(fp_psf), os.path.basename(fp_psfto))
	# hdus.writeto(fp_img_out, overwrite=overwrite)

	fits.PrimaryHDU(psf_out).writeto(fp_psf_out, overwrite=overwrite)
	if towrite_psk:
		fits.PrimaryHDU(psk_out).writeto(fp_psk_out, overwrite=overwrite)


def replace_img_in_fits(fn_from, fn_to, img, comment='', overwrite=True):
	"""
	replace the image in fn_from fits file by img and save it to fn+to
	"""
	hdus = fits.open(fn_from)
	hdus[0].data = img
	if len(comment) > 0:
		hdus[0].header['COMMENT'] = comment
	hdus.writeto(fn_to, overwrite=overwrite)


def match_psf(img, psf, psfto): 
	"""
	match the psf of the image from psf to psfto

	Params
	------
	img (np 2d array)
	psf (np 2d array)
	psfto (np 2d array)

	Return
	------
	img_cnvled (np 2d array):
		the psf matched image
	psf_cnvled (np 2d array):
		the matched psf (psf convolved with kernel_cnvl)
	kernel_cnvl (np 2d array):
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

		p0 = pad_edge_to_shape(p0, int(nxmax), int(nymax))
		p1 = pad_edge_to_shape(p1, int(nxmax), int(nymax))

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
		arr = np.pad(arr, ((int(dx), int(dx)), (0, 0)), 'constant', constant_values=0)
	elif nx < nxa:
		raise Exception("[pad_edge_to_shape] final dimension smaller than array")

	if ny > nya:
		dy = (ny - nya)/2
		arr = np.pad(arr, ((0, 0), (int(dy), int(dy))), 'constant', constant_values=0)
	elif ny < nya:
		raise Exception("[pad_edge_to_shape] final dimension smaller than array")

	return arr


def calc_max_frac_diff_between_two_psf(fp_psf1, fp_psf2):
	"""
	calculate the maximum difference between the two psfs.
	The difference is shown as a fraction of the maximum of the maximum of the two normalized psfs.
	"""
	nx, ny = 43, 43
	psf1 = fits.getdata(fp_psf1) 
	psf2 = fits.getdata(fp_psf2)

	psf1 = normalize_kernel(psf1)
	psf2 = normalize_kernel(psf2)

	psf1 = pad_edge_to_shape(psf1, nx=nx, ny=ny)
	psf2 = pad_edge_to_shape(psf2, nx=nx, ny=ny)

	max_orig = max([psf1.max(), psf2.max()])
	max_diff = np.max(np.absolute(psf1-psf2))

	frac_diff = max_diff/max_orig
	return frac_diff


def normalize_kernel(kernel):
	return kernel/np.sum(kernel)


#======================== about fitting model to psf and psf size measure =================

def get_xy_grid(nx, ny):
	"""
	return x_grid, y_grid symmetric and centered on central pixel
	accepting only odd nx, ny
	"""
	for n in [nx, ny]:
		if not isodd(n):
			raise Exception("[get_xy_grid] only accept odd number")

	x, y = np.mgrid[-(nx-1)/2:(nx+1)/2, -(ny-1)/2:(ny+1)/2]

	return x, y



def fit_moffat(arr):
	"""
	Params
	------
	arr (2d np array of odd sizes)

	Return
	------
	model (astropy model object)
	"""

	# sanity check of input data type
	if isinstance(arr, ac.kernels.Kernel):
		arr = arr.array
	elif isinstance(arr, np.ndarray):
		pass
	else: 
		raise Exception("[psfmatch] input needs to be a kernel or array")

	nx, ny = arr.shape
	x, y = get_xy_grid(nx, ny)

	model_init = am.functional_models.Moffat2D(amplitude=arr.max(), x_0=0, y_0=0, gamma=1, alpha=1)
	fitter = am.fitting.LevMarLSQFitter()

	# with warnings.catch_warnings():
		# warnings.simplefilter('ignore')
	model_best = fitter(model_init, x, y, arr)

	if model_best.gamma < 0:
		model_best.gamma = -model_best.gamma

	return model_best


def fit_gaussian(arr):
	"""
	Params
	------
	arr (2d np array of odd sizes)

	Return
	------
	model (astropy model object)
	"""
	if isinstance(arr, ac.kernels.Kernel):
		arr = arr.array
	elif isinstance(arr, np.ndarray):
		pass
	else: 
		raise Exception("[psfmatch] input needs to be a kernel or array")

	nx, ny = arr.shape
	x, y = get_xy_grid(nx, ny)

	model_init = am.functional_models.Gaussian2D(amplitude=arr.max(), x_mean=0., y_mean=0., x_stddev=5., y_stddev=5., theta=0.)
	fitter = am.fitting.LevMarLSQFitter()

	# with warnings.catch_warnings():
		# warnings.simplefilter('ignore')
	model_best = fitter(model_init, x, y, arr)

	return model_best



def calc_psf_fwhm(arr, mode='moffat'):
	""" 
	calculate the fwhm of psf by fitting the psf with model (either moffat or gaussian) and then calculate the analytic fwhm based on the params of the best fit. fwhm is expressed in pixel units. 

	Params
	------
	arr (2d np arr)
	mode (str):
		moffat, gaussian

	Return
	------
	fwhm (float)
		fwhm of psf in pix
	"""

	if mode == 'moffat':
		return calc_psf_fwhm_inpix_moffat(arr)
	elif mode == 'gaussian':
		return calc_psf_fwhm_inpix_gaussian(arr)
	else:
		raise ValueError("mode not recognized")


def calc_psf_fwhm_inpix_moffat(arr):
	""" 
	Params
	------
	arr (2d np arr)

	Return
	------
	fwhm (float)
		fwhm of psf in pix
	"""
	model = fit_moffat(arr)

	fwhm = 2.* model.gamma * np.sqrt( 2.**(1./model.alpha) - 1. )

	return fwhm


def calc_psf_fwhm_inpix_gaussian(arr):
	""" 
	takes the fwhm of the major axis

	Params
	------
	arr (2d np arr)

	Return
	------
	fwhm (float)
		fwhm of psf in pix
	"""
	model = fit_gaussian(arr)

	sigma = max(model.y_stddev, model.x_stddev)
	fwhm = 2.355 * sigma

	return fwhm


def calc_psf_size_inpix_quick(arr):
	""" quick and dirty way to calculate the size of the psf with arbitrary scaling, just for comparison. about 5 times faster than the other methods """
	arr1d = arr.sum(axis=0)
	x = np.arange(arr1d.size)
	spline = si.UnivariateSpline(x, arr1d-np.max(arr1d)/2, s=0)
	r1, r2 = spline.roots()

	return np.absolute(r2 - r1)


def has_smaller_psf_fits(fp1, fp2, mode='quick', fracdiff_threshold=0.1):
	""" 
	whether psf1 is smaller than psf2

	Params
	------
	fp1 (str):
		filepath to psf fits
	fp2 (str):
		filepath to psf fits
	mode (str):
		quick, moffat, or gaussian
	fracdiff_threshold (float):
		the fractional difference of psf sizes only beyond which will it be qualified as 'smaller'. 

	Return
	------
	results (bool)
	"""

	psf1 = fits.getdata(fp1)
	psf2 = fits.getdata(fp2)

	return has_smaller_psf(psf1, psf2, mode=mode, fracdiff_threshold=fracdiff_threshold)


def has_smaller_psf(psf1, psf2, mode='quick', fracdiff_threshold=0.1):
	""" given two psf 2d arrays determine if psf1 is smaller than psf2 """
	if mode == 'quick': 
		size1 = calc_psf_size_inpix_quick(psf1)
		size2 = calc_psf_size_inpix_quick(psf2)
	elif mode in ['moffat', 'gaussian']:
		size1 = calc_psf_fwhm(psf1, mode=mode)
		size2 = calc_psf_fwhm(psf2, mode=mode)
	else:
		raise ValueError("mode not recognized")


	return (size2 - size1)/size1 > fracdiff_threshold
	# return (size1 < size2)
