# test_psf.py
# ALS 2017/05/08
# testing set for psf.py


import pytest
import astropy.table as at
import os
import shutil

from ....obsobj import obsObj
from .. import psf


ra = 150.0547735
dec = 12.7073027

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent = './testing/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after test"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


@pytest.fixture
def obj_j1000():
	obj = obsObj(ra=ra, dec=dec, dir_parent=dir_parent)
	obj.add_sdss()

	return obj


def test_download_psField(obj_j1000):
	obj = obj_j1000

	psf.download_psField(xid=obj.sdss.xid, dir_out=obj.dir_obj)

	assert os.path.isfile(obj.dir_obj+'psField.fits')


def test_psField_to_psfs(obj_j1000):
	obj = obj_j1000

	psf.download_psField(xid=obj.sdss.xid, dir_out=obj.dir_obj)

	assert os.path.isfile(obj.dir_obj+'psField.fits')

	psf.psField_to_psfs(dir_obj=obj.dir_obj, photoobj=obj.sdss.photoobj, bands=['u', 'g', 'r', 'i', 'z'])

	for band in ['u', 'g', 'r', 'i', 'z']:
		assert os.path.isfile(obj.dir_obj+'psf-{0}.fits'.format(band))
