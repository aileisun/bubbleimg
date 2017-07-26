
import pytest
import os
import shutil
import numpy as np
import skimage.io as si

from .. import visualtools

dir_test = './testing/'
dir_verif = 'verification_data/'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_test):
		shutil.rmtree(dir_test)

	shutil.copytree(dir_verif, dir_test)

	yield
	# tear down
	if os.path.isdir(dir_test):
		shutil.rmtree(dir_test)


def test_visualtools_fits_to_png():

	fn_in = dir_test+'stamp-OIII5008_I.fits'
	visualtools.fits_to_png(fn_in, scaling='linear')

	assert os.path.isfile(dir_test+'stamp-OIII5008_I.png')


	fn_out = dir_test+'stamp-OIII5008_I_arcsinh.png'
	visualtools.fits_to_png(fn_in, fn_out=fn_out, scaling='arcsinh')
	assert os.path.isfile(fn_out)


	fn_out = dir_test+'stamp-OIII5008_I_arcsinh_trim.png'
	visualtools.fits_to_png(fn_in, fn_out=fn_out, scaling='arcsinh', vmax=10.)
	assert os.path.isfile(fn_out)



	fn_out = dir_test+'stamp-OIII5008_I_trim.png'
	visualtools.fits_to_png(fn_in, fn_out=fn_out, scaling='linear', vmax=10.)
	assert os.path.isfile(fn_out)



def test_visualtoools_scale_img():

	epsilon = 1.e-3

	img = np.zeros([3, 3])
	img[0, 0] = 1.

	img_scaled = visualtools.scale_img(img, vmin=None, vmax=None)

	img[0, 0] = 1. - epsilon

	assert img_are_similar(img_scaled, img[::-1, :])


	img2 = np.zeros([3, 3])
	img2[0, 0] = 2.

	img_scaled = visualtools.scale_img(img2, vmin=None, vmax=None)

	assert img_are_similar(img_scaled, img[::-1, :])

	fn_out = dir_test+'testing.png'
	si.imsave(fn_out, img_scaled)



def img_are_similar(img1, img2, threshold = 1.e-3):

	return np.all(np.absolute(img1 - img2) < threshold)