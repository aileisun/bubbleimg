import os
import shutil
import numpy as np
import pytest

from astropy.io import fits
import astropy.convolution as ac

from .. import matchpsf


b0 = 'i'
b1 = 'z'
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'
dir_verif = 'test_verification_data/SDSSJ0920+0034/'


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


def test_deconvl_psf_data():

	psf0 = fits.getdata(dir_obj+'psf-{}.fits'.format(b0))
	psf1 = fits.getdata(dir_obj+'psf-{}.fits'.format(b1))

	kernel = matchpsf.get_convolving_moffat_kernel(psf0, psf1)

	psf1_trial = ac.convolve(psf0, kernel)

	assert two_kernels_are_similar(psf1, psf1_trial)


def test_deconvl_match_psf():
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


