# test_loader_sdss.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for loader_sdss
"""

import numpy as np
import astropy.units as u
import shutil
import os
import pytest
from astropy.io import fits
import filecmp

from ..loader_sdss import SDSSimgLoader
# from ...class_obsobj import obsobj

ra = 150.0547735
dec = 12.7073027
img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent = './testing/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after test"""

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


@pytest.fixture
def L_radec_64pix():
	""" returns a imgLoader object initiated with the ra dec above of img size 64*64"""
	return SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64, img_height=64)


def test_make_stamps_write_frame(L_radec):
	""" test that frame can be downloaded and kept if tokeepframe=True, and deleted if tokeepframe=False """
	
	L = L_radec

	L.make_stamps(overwrite=True, tokeepframe=True)

	for band in L.bands:
		assert os.path.isfile(dir_obj+'frame-{0}.fits'.format(band))

	L.make_stamps(overwrite=True, tokeepframe=False)

	for band in L.bands:
		assert not os.path.isfile(dir_obj+'frame-{0}.fits'.format(band))


def test_make_stamps_overwriteTrue(L_radec):
	"""
	test that when overwrite=True make_stamp() always call download_stamp() whether file exists or not
	"""
	L = L_radec

	band = 'r'
	overwrite = True

	file = dir_obj+'stamp-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it creates stamp
	assert not os.path.isfile(file)
	L.make_stamps(overwrite=overwrite)
	assert os.path.isfile(file)

	# when file does exist it overwrites stamp
	if os.path.isfile(file):
		os.remove(file)
	open(file, 'w').close()

	L.make_stamps(overwrite=overwrite)
	assert os.path.isfile(file)
	assert os.stat(file).st_size > 0


def test_make_a_zero_size_file():
	file = 'testfile.txt'
	open(file, 'w').close()
	assert os.stat(file).st_size == 0
	os.remove(file)


def test_make_stamps_overwriteFalse(L_radec):
	"""
	test that when overwrite=False make_stamp() does not update file
	"""
	L = L_radec

	overwrite = False

 	for band in L.bands:
		file = dir_obj+'stamp-{0}.fits'.format(band)

		if os.path.isfile(file):
			os.remove(file)
		open(file, 'w').close()
		assert os.stat(file).st_size == 0

	# when file exists it should not update file
	L.make_stamps(overwrite=overwrite)

 	for band in L.bands:
		file = dir_obj+'stamp-{0}.fits'.format(band)
		assert os.stat(file).st_size == 0


def test_make_stamps_correctimgsize(L_radec_64pix):
	L = L_radec_64pix
	L.make_stamps(tokeepframe=False, overwrite=True)
	data = fits.getdata(L.dir_obj+'stamp-r.fits')
	assert data.shape == (64, 64)


def test_make_stamps_correctcontent(L_radec_64pix):
	L = L_radec_64pix
	L.make_stamps(tokeepframe=False, overwrite=True)

 	for band in L.bands:
 		f = 'stamp-{0}.fits'.format(band)
		file_totest = dir_obj+f
		file_verification = './test_verification_data_64pix/SDSSJ1000+1242/'+f
		assert filecmp.cmp(file_totest, file_verification)

	for f in ['sdss_xid.csv', 'sdss_photoobj.csv']:
		file_totest = dir_obj+f
		file_verification = './test_verification_data_64pix/SDSSJ1000+1242/'+f
		assert filecmp.cmp(file_totest, file_verification)


def test_get_stamp(L_radec_64pix):
	L = L_radec_64pix

	data = L.get_stamp(band='r')

	assert data.shape == (64, 64) 
	assert np.sum(data) > 0 
