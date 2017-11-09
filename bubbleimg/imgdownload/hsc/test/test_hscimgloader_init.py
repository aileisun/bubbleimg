# test_hscimgloader_init.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for hscimgloader

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

from ..hscimgloader import hscimgLoader
from ....obsobj import obsObj


ra     = 140.099364238908123 
dec    = 0.580160150759375104

img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './testing/SDSSJ0920+0034/'
dir_parent1 = './testing/'
dir_parent2 = './testing2/'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./testing2/ before and after test"""

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
	""" returns a hscimgLoader object initiated with the ra dec above"""
	return hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)


def test_automatically_generate_dir_obj_w_SDSSNAME():
	L = hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent1, img_width=img_width, img_height=img_height)
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_radec(L_radec):
	"""
	test that hscimgLoader can be instantiated with ra dec 
	"""

	L = L_radec

	assert isinstance(L, hscimgLoader)
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_obsobj():
	"""
	test that hscimgLoader can be instantiated with obsobj
	"""
	
	obj = obsObj(ra=ra, dec=dec, dir_parent=dir_parent2)

	L = hscimgLoader(obj=obj, img_width=img_width, img_height=img_height)

	assert isinstance(L, hscimgLoader)
	assert L.ra == ra
	assert L.ra == obj.ra
	assert L.img_height == img_height
	assert L.dir_obj == obj.dir_obj
	assert L.ra == L.obj.ra
	assert L.dec == L.obj.dec


def test_instantiate_HSCimgLoader_error_radec_obsobj():
	"""
	test that an error being raised when both obsobj and ra/dec/dir_obj are fed to hscimgLoader
	"""
	
	obj = obsObj(ra=ra, dec=dec, dir_parent=dir_parent2)

	with pytest.raises(Exception):
		L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj+'fortesting/', obj=obj, img_width=img_width, img_height=img_height)


def test_instantiate_HSCimgLoader_error():
	"""
	test that an error being raised when none of obsobj or ra/dec/dir_obj are fed to hscimgLoader
	"""

	with pytest.raises(TypeError):
		L = hscimgLoader(img_width=img_width, img_height=img_height)


def test_init_survey(L_radec):

	L = L_radec

	assert L.survey == 'hsc'
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_add_obj_sdss(L_radec):
	"""
	test that the function add_obj_sdss adds/updates L.obj property based on L.ra, dec, dir_obj, and also properly downloads xid.csv and photoobj.csv
	"""
	L = L_radec

	status = L.add_obj_sdss(update=False)
	assert status
	assert hasattr(L, 'obj')
	assert round(L.obj.ra, 3) == round(L.ra, 3)
	assert round(L.obj.dec, 3) == round(L.dec, 3)
	assert L.obj.dir_obj == L.dir_obj
	assert hasattr(L.obj, 'sdss')
	assert hasattr(L.obj.sdss, 'xid')
	assert hasattr(L.obj.sdss, 'camcol')
	assert os.path.isfile(L.dir_obj+'sdss_xid.csv')
	assert os.path.isfile(L.dir_obj+'sdss_photoobj.csv')
	xid = L.obj.sdss.xid

	# check that L.obj can be updated to the right things
	L.obj = 'testing'
	L.add_obj_sdss(update=True)

	assert round(L.obj.ra, 3) == round(L.ra, 3)
	assert L.obj.dir_obj == L.dir_obj
	assert L.obj.sdss.xid == xid



def test_add_obj_hsc(L_radec):
	"""
	test that the function add_obj_hsc adds/updates L.obj property and also properly downloads hsc_xid.csv
	"""
	L = L_radec

	assert L.status
	assert hasattr(L, 'obj')
	assert round(L.obj.ra, 3) == round(L.ra, 3)
	assert round(L.obj.dec, 3) == round(L.dec, 3)
	assert L.obj.dir_obj == L.dir_obj
	assert hasattr(L.obj, 'hsc')
	assert hasattr(L.obj.hsc, 'xid')
	assert hasattr(L.obj.hsc, 'tract')
	assert os.path.isfile(L.dir_obj+'hsc_xid.csv')
	assert L.obj.hsc.tract == 9564

	xid = L.obj.hsc.xid

	# check that L.obj can be updated to the right things
	L.obj = 'testing'
	L.add_obj_hsc(update=True)

	assert round(L.obj.ra, 3) == round(L.ra, 3)
	assert L.obj.dir_obj == L.dir_obj
	assert L.obj.hsc.xid == xid


def test_instantiate_HSCimgLoader_floatwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64., img_height=64)

	assert isinstance(L, hscimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_instantiate_HSCimgLoader_pixwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64.*u.pix, img_height=64.*u.pix)

	assert isinstance(L, hscimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_transform_img_widthheight_unit_to_pix(L_radec):
	L = L_radec

	hscpixsize = 0.168*u.arcsec/u.pix

	assert L.img_width_pix == int(np.floor((L.img_width/hscpixsize).to(u.pix).value))
	assert L.img_height_pix == int(np.floor((L.img_height/hscpixsize).to(u.pix).value))


def test_HSCimgLoader_get_img_width_pix(L_radec): 
	L = L_radec
	print(L.img_width_pix)
	assert type(L.img_width_pix) is int
	assert (L.img_width_pix) == 119.


def test_HSCimgLoader_get_img_width_arcsec(): 
	
	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=120, img_height=120)
	
	assert round(L.img_width_arcsec.value, 3) == 20.160
	assert round(L.img_height_arcsec.value, 3) == 20.160

	
