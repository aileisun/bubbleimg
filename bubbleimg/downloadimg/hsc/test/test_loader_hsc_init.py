# test_loader_hsc_init.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for loader_hsc

test suite init

please create a python file STARs.py with your STARs account
user = 'XXXX'
password = 'XXXX'

"""

import numpy as np
import astropy.table as at
import astropy.units as u
import shutil
import os

import pytest

from loader_hsc import HSCimgLoader
from ...class_obsobj import obsobj
from STARs import user, password

ra = 150.0547735
dec = 12.7073027
img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './test/SDSSJ1000+1242/'
dir_parent1 = './test/'
dir_parent2 = './test2/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./test/ and ./test2/ before and after test"""

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
	""" returns a HSCimgLoader object initiated with the ra dec above"""
	return HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height, user=user, password=password)


def test_automatically_generate_dir_obj_w_SDSSNAME():
	L = HSCimgLoader(ra=ra , dec=dec, dir_parent=dir_parent1, img_width=img_width, img_height=img_height, user=user, password=password)
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_radec(L_radec):
	"""
	test that HSCimgLoader can be instantiated with ra dec 
	"""

	L = L_radec

	assert isinstance(L, HSCimgLoader)
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_obsobj():
	"""
	test that HSCimgLoader can be instantiated with obsobj
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent=dir_parent2, towriteID=False)

	L = HSCimgLoader(obj=obj, img_width=img_width, img_height=img_height, user=user, password=password)

	assert isinstance(L, HSCimgLoader)
	assert L.ra == ra
	assert L.ra == obj.ra
	assert L.img_height == img_height
	assert L.dir_obj == obj.dir_obj
	assert L.ra == L.obj.ra
	assert L.dec == L.obj.dec


def test_instantiate_HSCimgLoader_error_radec_obsobj():
	"""
	test that an error being raised when both obsobj and ra/dec/dir_obj are fed to HSCimgLoader
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent=dir_parent2, towriteID=False)

	with pytest.raises(TypeError):
		L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, obj=obj, img_width=img_width, img_height=img_height, user=user, password=password)


def test_instantiate_HSCimgLoader_error():
	"""
	test that an error being raised when none of obsobj or ra/dec/dir_obj are fed to HSCimgLoader
	"""

	with pytest.raises(TypeError):
		L = HSCimgLoader(img_width=img_width, img_height=img_height, user=user, password=password)


def test_init_survey(L_radec):

	L = L_radec

	assert L.survey == 'hsc'
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_add_obj_sdss(L_radec):
	"""
	test that the function _add_obj_sdss adds/updates L.obj property based on L.ra, dec, dir_obj, and also properly downloads xid.csv and photoobj.csv
	"""
	L = L_radec
	assert not hasattr(L, 'obj')

	L._add_obj_sdss(update=False)
	assert hasattr(L, 'obj')
	assert L.obj.ra == L.ra
	assert L.obj.dec == L.dec
	assert L.obj.dir_obj == L.dir_obj
	assert hasattr(L.obj, 'sdss')
	assert hasattr(L.obj.sdss, 'xid')
	assert os.path.isfile(L.dir_obj+'xid.csv')
	assert os.path.isfile(L.dir_obj+'photoobj.csv')
	xid = L.obj.sdss.xid

	# check that L.obj can be updated to the right things
	L.obj = 'testing'
	L._add_obj_sdss(update=True)

	assert L.obj.ra == L.ra
	assert L.obj.dir_obj == L.dir_obj
	assert L.obj.sdss.xid == xid


def test_make_obj_sdss_xid(L_radec):
	L = L_radec
	assert not hasattr(L, 'obj')

	L._make_obj_sdss_xid()
	assert round(L.obj.sdss.xid['ra'][0], 4) == round(L.ra, 4)
	assert os.path.isfile(L.dir_obj+'xid.csv')


def test_instantiate_HSCimgLoader_floatwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64., img_height=64, user=user, password=password)

	assert isinstance(L, HSCimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_pixwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64.*u.pix, img_height=64.*u.pix, user=user, password=password)

	assert isinstance(L, HSCimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_transform_img_widthheight_unit_to_pix(L_radec):
	L = L_radec

	hscpixsize = 0.168*u.arcsec/u.pix

	assert L.img_width_pix == np.floor((L.img_width/hscpixsize).to(u.pix))
	assert L.img_height_pix == np.floor((L.img_height/hscpixsize).to(u.pix))


def test_HSCimgLoader_get_img_width_pix(L_radec): 
	L = L_radec
	print L.img_width_pix
	assert (L.img_width_pix.value).is_integer()
	assert (L.img_width_pix.value) == 119.


def test_HSCimgLoader_get_img_width_arcsec(): 
	
	L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=120, img_height=120, user=user, password=password)
	
	assert round(L.img_width_arcsec.value, 3) == 20.160
	assert round(L.img_height_arcsec.value, 3) == 20.160

	
