# test_obsobj.py
# ALS 2017/05/16

import pytest
import os
import shutil

from ..obsobj import obsObj

ra = 150.0547735
dec = 12.7073027

dir_obj = './test/SDSSJ1000+1242/'
dir_parent = './test/'


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
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


def test_obsObj_init_dir_obj(obj_dirobj):

	obj = obj_dirobj

	assert isinstance(obj, obsObj)
	assert obj.dir_obj == dir_obj


def test_obsObj_init_dir_parent():

	obj = obsObj(ra=ra, dec=dec, dir_parent = dir_parent)

	assert isinstance(obj, obsObj)
	assert obj.dir_obj == dir_obj


def test_obsObj_add_sdss(obj_dirobj):
	obj = obj_dirobj

	status = obj.add_sdss()

	assert status
	assert round(obj.sdss.ra, 4) == round(obj.ra, 4)
	assert round(obj.sdss.dec, 4) == round(obj.dec, 4)
	assert hasattr(obj.sdss, 'xid')
	assert hasattr(obj.sdss, 'photoobj')



def test_obsObj_add_sdss_fail():

	ra = 0.
	dec = -89.
	obj = obsObj(ra=ra, dec=dec, dir_parent = dir_parent)

	status = obj.add_sdss()
	assert status == False

	assert obj.sdss.status == False
