import os
import shutil
import numpy as np
import pytest

import astropy.table as at
from astropy.io import fits
import astropy.convolution as ac

from .. import matchpsf


b0 = 'i'
b1 = 'z'
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'
dir_verif = 'verification_data/SDSSJ0920+0034/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	os.makedirs(dir_parent)
	shutil.copytree(dir_verif, dir_obj)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


def two_kernels_are_similar(k0, k1, avgdiff_threshold=1.e-3):

	try:
		k0 = k0.array
	except:
		pass

	try:
		k1 = k1.array
	except:
		pass

	k0_m, k1_m = matchpsf.match_dimension(k0, k1)

	avgdiffsq = np.average((k0_m - k1_m)**2, weights=k0_m)

	avgdiff = np.sqrt(avgdiffsq)
	print("avgdiff "+str(avgdiff))
	return avgdiff < avgdiff_threshold


def test_kernels_are_similar():

	psfi = fits.getdata(dir_obj+'psf-i.fits')
	psfr = fits.getdata(dir_obj+'psf-r.fits')
	
	# avgdiff = 0.0180033021648
	assert not two_kernels_are_similar(psfi, psfr)


def test_pad_edge_to_shape():
	nx, ny = 5, 7
	arr = np.ndarray([1, 3])

	arr_pad = matchpsf.pad_edge_to_shape(arr, nx, ny)

	assert arr_pad.shape == (nx, ny)
	assert arr.shape == (1, 3)
	assert arr_pad[0, 0] == 0. 

	with pytest.raises(Exception) as e:
		nx, ny = 5, 6
		arr_pad = matchpsf.pad_edge_to_shape(arr, nx, ny)

	with pytest.raises(Exception) as e:
		nx, ny = 5, 7
		arr = np.ndarray([2, 3])
		arr_pad = matchpsf.pad_edge_to_shape(arr, nx, ny)


def test_match_dimension():

	p0 = np.ndarray([1, 3])
	p1 = np.ndarray([5, 7])
	p0, p1 = matchpsf.match_dimension(p0, p1)

	assert p0.shape == p1.shape

	with pytest.raises(Exception) as e:
		p0 = np.ndarray([1, 3])
		p1 = np.ndarray([6, 7])
		p0, p1 = matchpsf.match_dimension(p0, p1)


def test_is_odd():
	assert matchpsf.isodd(1)
	assert not matchpsf.isodd(2)

def test_normalize_kernel():
	kernel = np.random.rand(3, 3)

	kernel_norm = matchpsf.normalize_kernel(kernel)
	assert round(np.sum(kernel_norm), 7) == 1.


def test_psf_data():

	psf0 = fits.getdata(dir_obj+'psf-{}.fits'.format(b0))
	psf1 = fits.getdata(dir_obj+'psf-{}.fits'.format(b1))

	kernel = matchpsf.get_convolving_moffat_kernel(psf0, psf1)

	psf1_trial = ac.convolve(psf0, kernel)

	assert two_kernels_are_similar(psf1, psf1_trial)


def test_match_psf():
	""" test that match_psf can infer the simulated moffat kernel """

	band = b0

	img = fits.getdata(dir_obj+'stamp-{}.fits'.format(band))
	psf = fits.getdata(dir_obj+'psf-{}.fits'.format(band))

	kernel = ac.Moffat2DKernel(gamma=5., alpha=3., x_size=43, y_size=43)

	psfto = ac.convolve(psf, kernel)
	img_cnvled_veri = ac.convolve(img, kernel)

	img_cnvled, psf_cnvled, kernel_cnvl = matchpsf.match_psf(img, psf, psfto)


	assert two_kernels_are_similar(kernel_cnvl, kernel)
	assert two_kernels_are_similar(psf_cnvled, psfto)
	assert two_kernels_are_similar(img_cnvled, img_cnvled_veri)




def test_match_psf_fits():

	bandfrom = b0
	bandto = b1

	fp_img = dir_obj+'stamp-{}.fits'.format(bandfrom)
	fp_psf = dir_obj+'psf-{}.fits'.format(bandfrom)
	fp_psfto = dir_obj+'psf-{}.fits'.format(bandto)

	fp_img_out = dir_obj+'stamp-{}_psfmt-{}.fits'.format(bandfrom, bandto)
	fp_psf_out = dir_obj+'psf-{}_psfmt-{}.fits'.format(bandfrom, bandto)
	fp_psk_out = dir_obj+'psk-{}_psfmt-{}.fits'.format(bandfrom, bandto)

	matchpsf.match_psf_fits(fp_img, fp_psf, fp_psfto, fp_img_out, fp_psf_out, fp_psk_out, overwrite=True, towrite_psk=True)	

	# check that file exists
	for fp in [fp_img_out,fp_psf_out,fp_psk_out,]:
		assert os.path.isfile(fp)


	img = fits.getdata(fp_img)
	psf = fits.getdata(fp_psf)
	psfto = fits.getdata(fp_psfto)

	# check that the final image content is correct
	img_out, psf_out, psk_out = matchpsf.match_psf(img, psf, psfto)

	assert np.all(img_out == fits.getdata(fp_img_out))
	assert np.all(psf_out == fits.getdata(fp_psf_out))



def test_get_xy_grid():
	nx = 3
	ny = 5
	x, y = matchpsf.get_xy_grid(nx, ny)

	assert x[0,0] == - x[-1,0]
	assert y[0,0] == - y[0,-1]
	assert x.shape == (nx, ny)
	assert y.shape == (nx, ny)

	with pytest.raises(Exception) as e:
		x, y = matchpsf.get_xy_grid(4, 6)


@pytest.mark.parametrize("fit_psf_func, frac_peak_err", [(matchpsf.fit_moffat, 0.05), (matchpsf.fit_gaussian, 0.5)])
def test_fit_psf_func(fit_psf_func, frac_peak_err):
	"""fit psf with moffat and return best fit model
	frac_peak_err is the tolerable error as a fraction of the peak psf value 
	"""

	psf = fits.getdata(dir_verif+'psf-r.fits')
	x, y = matchpsf.get_xy_grid(*psf.shape)

	model = fit_psf_func(psf)

	psf_model = model(x, y)

	diff = (psf-psf_model)/psf
	diff[~np.isfinite(diff)] = 0.

	assert np.absolute(psf-psf_model).max() < frac_peak_err*psf.max()


def test_calc_psf_fwhm_inpix():
	psf = fits.getdata(dir_verif+'psf-r.fits')


	for mode in ['moffat', 'gaussian']:
		fwhm_inpix = matchpsf.calc_psf_fwhm(psf, mode=mode)

		# fwhm_inpix = matchpsf.calc_psf_fwhm_inpix_moffat(psf)

		fwhm_manual_measured_inpix = 4.869

		fracerr = np.absolute((fwhm_inpix - fwhm_manual_measured_inpix)/fwhm_manual_measured_inpix)
		print fracerr
		assert fracerr < 0.2

	# moffat is more accurate




def test_has_smaller_psf_fits():

	for mode in ['quick', 'moffat', 'gaussian']:
		assert matchpsf.has_smaller_psf_fits(dir_verif+'psf-i.fits', dir_verif+'psf-r.fits', mode=mode)
		assert matchpsf.has_smaller_psf_fits(dir_verif+'psf-i.fits', dir_verif+'psf-y.fits', mode=mode)

