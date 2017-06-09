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

from ..hscimgloader import hscimgLoader

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


@pytest.fixture(params=['iaa', 'online'])
def L_fixture(request):
	""" returns a imgLoader object initiated with the ra dec above"""
	return hscimgLoader(ra=ra, dec=dec, dir_parent=dir_parent, environment=request.param)


def test_make_psf_filename(L_fixture):
	L = L_fixture
	fn = L.get_fn_psf(band='z')
	assert fn == 'psf-z.fits'


def test_make_psf(L_fixture):
	""" test it can make_stamp"""
	L = L_fixture
	
	status = L.make_psf(band = 'r', overwrite=True)

	assert status
	fn = L.dir_obj+'psf-r.fits'
	assert os.path.isfile(fn)

	fn_veri = './test_verification_data_128pix/SDSSJ0920+0034/psf-r.fits'

	# check that psf is the same as verfication
	data = fits.getdata(fn)
	data_veri = fits.getdata(fn_veri)
	assert np.max(np.absolute(data-data_veri)) < 0.0001
	# assert filecmp.cmp(fn, fn_veri)


def test_make_psf_not_keep_calexp(L_fixture):

	L = L_fixture
	status = L.make_psf(band = 'r', overwrite=True, to_keep_calexp=False)

	assert status
	assert os.path.isfile(L.dir_obj+'psf-r.fits')
	assert not os.path.isfile(L.dir_obj+'calexp-r.fits')


def test_make_psf_distinct_psf(L_fixture):

	L = L_fixture
	status = L.make_psf(band = 'r', overwrite=True, to_keep_calexp=False)
	assert status
	status = L.make_psf(band = 'z', overwrite=True, to_keep_calexp=False)
	assert status

	# the two files are distinct
	assert not filecmp.cmp(L.dir_obj+'psf-r.fits', L.dir_obj+'psf-z.fits')


def test_make_psf_overwrite(L_fixture):
	# when the file exists it still overwrites
	L = L_fixture

	fn = L.dir_obj+'psf-r.fits'

	if os.path.isfile(fn):
		os.remove(fn)
	open(fn, 'w').close()

	assert os.stat(fn).st_size == 0

	status = L.make_psf(band = 'r', overwrite=True, to_keep_calexp=False)

	assert os.path.isfile(fn)
	assert os.stat(fn).st_size > 0


def test_make_psfs_overwrite(L_fixture):
	# when the file exists it still overwrites
	L = L_fixture

	for band in L.bands:
		fn = L.dir_obj+'psf-{}.fits'.format(band)

		if os.path.isfile(fn):
			os.remove(fn)
		open(fn, 'w').close()

		assert os.stat(fn).st_size == 0

	status = L.make_psfs(overwrite=True, to_keep_calexp=False)

	for band in L.bands:
		fn = L.dir_obj+'psf-{}.fits'.format(band)
		assert os.path.isfile(fn)
		assert os.stat(fn).st_size > 0


@pytest.mark.parametrize("environment", [('iaa'), ('online')])
def test_make_psf_no_img_in_hsc(environment):
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, environment=environment)
	status = L.make_psf(band = 'r', overwrite=True)
	assert status is False
