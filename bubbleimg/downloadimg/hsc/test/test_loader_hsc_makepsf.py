# test_loader_makepsf.py
# ALS 2017/05/23

import numpy as np
import astropy.units as u
import shutil
import os
import pytest
from astropy.io import fits
import filecmp
# import glob

from ..loader_hsc import HSCimgLoader

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'


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
	return HSCimgLoader(ra=ra, dec=dec, dir_parent=dir_parent)


def test_make_psf(L_radec):
	""" test it can make_stamp"""
	L = L_radec
	status = L.make_psf(band = 'r', overwrite=True)

	fn = L.dir_obj+'psf-r.fits'
	assert status
	assert os.path.isfile(fn)

	file_verification = './test_verification_data_128pix/SDSSJ0920+0034/psf-r.fits'
	assert filecmp.cmp(fn, file_verification)


def test_make_psf_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj)

	status = L.make_psf(band = 'r', overwrite=True)
	assert status is False


def test_make_psf_not_keep_calexp(L_radec):

	L = L_radec
	status = L.make_psf(band = 'r', overwrite=True, to_keep_calexp=False)

	assert status
	assert os.path.isfile(L.dir_obj+'psf-r.fits')
	assert not os.path.isfile(L.dir_obj+'calexp-r.fits')
