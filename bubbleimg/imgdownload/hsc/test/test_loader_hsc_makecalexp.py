# test_loader_hsc_makecalexp.py
# ALS 2017/05/24

"""
to be used with pytest

test sets for loader_hsc
"""

import numpy as np
import astropy.units as u
import shutil
import os
import pytest
from astropy.io import fits
import filecmp
import glob

from ..loader_hsc import hscimgLoader

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'
img_width = 128*u.pix
img_height = 128*u.pix


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
	return hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=img_width, img_height=img_height)


def test_make_calexp(L_radec):
	""" test it can make_calexp"""
	L = L_radec
	assert L.status
	
	fn = L.dir_obj+'calexp-r.fits'
	status = L.make_calexp(band = 'r', overwrite=True)

	assert status
	assert os.path.isfile(fn)

	file_verification = './test_verification_data_128pix/SDSSJ0920+0034/calexp-r.fits'
	assert filecmp.cmp(fn, file_verification)


def test_make_calexp_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	status = L.make_calexp(band = 'r', overwrite=True)
	assert status is False


def test_make_calexps(L_radec):
	""" test it can make_calexp"""
	L = L_radec
	status = L.make_calexps(overwrite=True)

	assert status

	for b in ['g', 'r', 'i', 'z', 'y']:
		assert os.path.isfile(L.dir_obj+'calexp-{0}.fits'.format(b))


