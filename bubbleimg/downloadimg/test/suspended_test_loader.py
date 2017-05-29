# test_imgLoader.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for imgLoader
"""

import numpy as np
import astropy.table as at
import astropy.units as u
import shutil
import os

import pytest

from ..loader import imgLoader

from ...class_obsobj import obsobj


ra = 150.0547735
dec = 12.7073027
img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent1 = './testing/'
dir_parent2 = './testing2/'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after test"""

	# setup
	if os.path.isdir(dir_parent1):
		shutil.rmtree(dir_parent1)

	if os.path.isdir(dir_parent2):
		shutil.rmtree(dir_parent2)

	yield
	# tear down
	if os.path.isdir(dir_parent1):
		shutil.rmtree(dir_parent1)

	if os.path.isdir(dir_parent2):
		shutil.rmtree(dir_parent2)


@pytest.fixture
def L_radec():
	""" returns a imgLoader object initiated with the ra dec above"""
	return imgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)



def test_make_stamp_mkdir(L_radec):
	"""
	test that make_stamp creates dir_obj when it's absent
	"""
	L = L_radec

	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)

	with pytest.raises(NotImplementedError):
		L.make_stamp(overwrite=True)

	assert os.path.isdir(dir_obj)
	shutil.rmtree(dir_obj)


def test_make_stamp_call_download_stamp(L_radec):
	"""
	test that when overwrite=True make_stamp() always call download_stamp() whether file exists or not
	"""
	L = L_radec

	band = 'r'
	overwrite = True

	file = dir_obj+'stamp-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it should call download_stamp
	with pytest.raises(NotImplementedError): 
		L.make_stamp(band=band, overwrite=overwrite)

	open(file, 'w').close()

	# when file exists it should call download_stamp
	with pytest.raises(NotImplementedError): 
		L.make_stamp(band=band, overwrite=overwrite)


def test_make_stamp_skipcall_download_stamp(L_radec):
	"""
	test that when overwrite=False make_stamp() does not call download_stamp() when file exists
	"""
	L = L_radec

	band = 'r'
	overwrite = False

	file = dir_obj+'stamp-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it should call download_stamp()
	with pytest.raises(NotImplementedError): 
		L.make_stamp(band=band, overwrite=overwrite)

	open(file, 'w').close()

	# when file exists nothing should happen
	L.make_stamp(band=band, overwrite=overwrite)


def test_make_stamps_error(L_radec):
	"""
	test that make_stamps fails because self.survey = 'to be overwritten'
	"""
	L = L_radec

	with pytest.raises(KeyError):
		L.make_stamps()


