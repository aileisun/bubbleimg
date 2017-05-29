# test_obsobj.py
# ALS 2017/05/16

import pytest
import os
import shutil

from ..obsobj import obsObj

ra = 150.0547735
dec = 12.7073027

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent = './testing/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

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


@pytest.fixture
def obj_hscobj():
	ra = 140.099341430207
	dec = 0.580162492432517
	dir_obj = './testing/SDSSJ0920+0034/'
	return obsObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_obsObj_init_dir_obj(obj_dirobj):

	obj = obj_dirobj

	assert isinstance(obj, obsObj)
	assert obj.dir_obj == dir_obj
	assert hasattr(obj, 'ra')
	assert hasattr(obj, 'dir_obj')
	assert hasattr(obj, 'name')

	assert obj.name == obj.dir_obj.split('/')[-2]



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


def test_obsObj_add_hsc(obj_hscobj):
	obj = obj_hscobj

	status = obj.add_hsc()

	assert status
	assert round(obj.hsc.ra, 2) == round(obj.ra, 2)
	assert round(obj.hsc.dec, 2) == round(obj.dec, 2)
	assert hasattr(obj.hsc, 'xid')
	assert hasattr(obj.hsc, 'tract')
	assert hasattr(obj.hsc, 'patch')


def test_obsObj_add_hsc_fail():

	ra = 0.
	dec = -89.
	obj = obsObj(ra=ra, dec=dec, dir_parent = dir_parent)

	status = obj.add_hsc()
	assert status == False

	assert obj.hsc.status == False

