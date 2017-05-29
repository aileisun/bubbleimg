# test_loader_sdss_init.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for loader_sdss

test suite init
"""

import numpy as np
import astropy.table as at
import astropy.units as u
import shutil
import os

import pytest

from ..loader_sdss import SDSSimgLoader

from ....class_obsobj import obsobj


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
	""" returns a SDSSimgLoader object initiated with the ra dec above"""
	return SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)


def test_instantiate_SDSSimgLoader_radec(L_radec):
	"""
	test that SDSSimgLoader can be instantiated with ra dec 
	"""

	L = L_radec

	assert isinstance(L, SDSSimgLoader)
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_instantiate_SDSSimgLoader_obsobj():
	"""
	test that SDSSimgLoader can be instantiated with obsobj
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent=dir_parent2, towriteID=False)

	L = SDSSimgLoader(obj=obj, img_width=img_width, img_height=img_height)

	assert isinstance(L, SDSSimgLoader)
	assert L.ra == ra
	assert L.ra == obj.ra
	assert L.img_height == img_height
	assert L.dir_obj == obj.dir_obj
	assert L.ra == L.obj.ra
	assert L.dec == L.obj.dec


def test_instantiate_SDSSimgLoader_error_radec_obsobj():
	"""
	test that an error being raised when both obsobj and ra/dec/dir_obj are fed to SDSSimgLoader
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent=dir_parent2, towriteID=False)

	with pytest.raises(Exception):
		L = SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, obj=obj, img_width=img_width, img_height=img_height)


def test_instantiate_SDSSimgLoader_error():
	"""
	test that an error being raised when none of obsobj or ra/dec/dir_obj are fed to SDSSimgLoader
	"""

	with pytest.raises(TypeError):
		L = SDSSimgLoader(img_width=img_width, img_height=img_height)


def test_init_survey(L_radec):

	L = L_radec

	assert L.survey == 'sdss'
	assert L.ra == ra
	assert L.img_width == img_width
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj



def test_add_obj_sdss(L_radec):
	"""
	test that the function _add_obj_sdss adds/updates L.obj property based on L.ra, dec, dir_obj, and also properly downloads xid.csv and photoobj.csv
	"""
	L = L_radec

	del L_radec.obj
	assert not hasattr(L, 'obj')

	L.add_obj_sdss(update=False)
	assert hasattr(L, 'obj')
	assert L.obj.ra == L.ra
	assert L.obj.dec == L.dec
	assert L.obj.dir_obj == L.dir_obj
	assert hasattr(L.obj, 'sdss')
	assert hasattr(L.obj.sdss, 'xid')
	assert os.path.isfile(L.dir_obj+'sdss_xid.csv')
	assert os.path.isfile(L.dir_obj+'sdss_photoobj.csv')
	xid = L.obj.sdss.xid

	# check that L.obj can be updated to the right things
	L.obj = 'testing'
	L.add_obj_sdss(update=True)

	assert L.obj.ra == L.ra
	assert L.obj.dir_obj == L.dir_obj
	assert L.obj.sdss.xid == xid


def test_instantiate_SDSSimgLoader_floatwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64., img_height=64)

	assert isinstance(L, SDSSimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_instantiate_SDSSimgLoader_pixwidth():
	"""
	test img_width can be input with float (with units assumed to be pix)
	"""

	L = SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=64.*u.pix, img_height=64.*u.pix)

	assert isinstance(L, SDSSimgLoader)
	assert L.ra == ra
	assert L.img_width == 64*u.pix
	assert L.img_height == 64*u.pix
	assert L.dir_obj == dir_obj


def test_transform_img_widthheight_unit_to_pix(L_radec):
	L = L_radec

	sdsspixsize = 0.396*u.arcsec/u.pix

	assert L.img_width_pix == np.floor((L.img_width/sdsspixsize).to(u.pix))
	assert L.img_height_pix == np.floor((L.img_height/sdsspixsize).to(u.pix))


def test_SDSSimgLoader_get_img_width_pix(L_radec): 
	L = L_radec
	print L.img_width_pix
	assert (L.img_width_pix.value).is_integer()
	assert (L.img_width_pix.value) == 50.
