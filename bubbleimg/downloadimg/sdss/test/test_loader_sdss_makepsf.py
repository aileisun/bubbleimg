# test_loader_sdss_makepsf.py
# ALS 2017/05/08

"""
to be used with pytest

test sets for loader_sdss make_psf
"""

import numpy as np
import astropy.units as u
import shutil
import os
import pytest
from astropy.io import fits
import filecmp


from ..loader_sdss import SDSSimgLoader
from ....class_obsobj import obsobj

ra = 150.0547735
dec = 12.7073027
img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent = './testing/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./test/ and ./test2/ before and after test"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def L_radec():
	""" returns a imgLoader object initiated with the ra dec above"""
	return SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)


def test_make_psfs_write_psf(L_radec):
	L = L_radec
	L.make_psfs(overwrite=True, tokeepfield=True)

	assert os.path.isfile(dir_obj+'psField.fits')
	for band in L.bands:
		assert os.path.isfile(dir_obj+'psf-{0}.fits'.format(band))


def test_make_psfs_overwriteTrue(L_radec):
	"""
	test that when overwrite=True make_psfs() always overwrite
	"""
	L = L_radec

	band = 'r'
	overwrite = True

	file = dir_obj+'psf-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it creates file
	assert not os.path.isfile(file)
	L.make_psfs(overwrite=overwrite)
	assert os.path.isfile(file)

	# when file does exist it overwrites file
	if os.path.isfile(file):
		os.remove(file)
	open(file, 'w').close()

	L.make_psfs(overwrite=overwrite)
	assert os.path.isfile(file)
	assert os.stat(file).st_size > 0


def test_make_psfs_overwriteFalse(L_radec):
	"""
	test that when overwrite=True make_psfs() does not overwrite
	"""
	L = L_radec

	overwrite = False

 	for band in L.bands:
		file = dir_obj+'psf-{0}.fits'.format(band)

		if os.path.isfile(file):
			os.remove(file)
		open(file, 'w').close()
		assert os.stat(file).st_size == 0

	# when file exists it should not update file
	L.make_psfs(overwrite=overwrite)

 	for band in L.bands:
		file = dir_obj+'psf-{0}.fits'.format(band)
		assert os.stat(file).st_size == 0


def test_make_psfs_not_tokeepField(L_radec):
	L = L_radec

	f_psfield = dir_obj+'psField.fits'

	if os.path.isfile(f_psfield):
		os.remove(f_psfield)
	assert not os.path.isfile(f_psfield)

	L.make_psfs(overwrite=True, tokeepfield=False)

	assert not os.path.isfile(f_psfield)


def test_make_psfs_correctcontent(L_radec):
	L = L_radec
	L.make_psfs(overwrite=True)

 	for band in L.bands:
 		f = 'psf-{0}.fits'.format(band)
		file_totest = dir_obj+f
		file_verification = './test_verification_data_64pix/SDSSJ1000+1242/'+f
		assert filecmp.cmp(file_totest, file_verification)


def test_psf_area_sums_to_one(L_radec):
	L = L_radec

 	for band in L.bands:
	
		psf = L.get_psf(band)
		assert len(psf) > 0
		assert round(np.sum(psf), 3) == 1.
