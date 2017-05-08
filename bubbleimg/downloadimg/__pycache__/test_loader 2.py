# test_imgLoader.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for imgLoader
"""

import numpy as np
import astropy.units as u
import astropy.table as at

import pytest

from loader import imgLoader

from ..class_obsobj import obsobj


ra = 150.0547735
dec = 12.7073027
img_width = 20*u.arcsec
img_height = 20*u.arcsec

dir_obj = './test/SDSSJ1000+1242/'



def test_instantiate_imgLoader_radec():
	"""
	test that imgLoader can be instantiated with ra dec 
	"""

	L = imgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	assert isinstance(L, imgLoader)
	assert L.ra == ra
	assert L.img_height == img_height
	assert L.dir_obj == dir_obj


def test_instantiate_imgLoader_obsobj():
	"""
	test that imgLoader can be instantiated with obsobj
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent='./test2/', towriteID=False)

	L = imgLoader(obsobj=obj, img_width=img_width, img_height=img_height)

	assert isinstance(L, imgLoader)
	assert L.ra == ra
	assert L.ra == obj.ra
	assert L.img_height == img_height
	assert L.dir_obj == obj.dir_obj


def test_instantiate_imgLoader_error_radec_obsobj():
	"""
	test that an error being raised when both obsobj and ra/dec/dir_obj are fed to imgLoader
	"""
	
	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent='./test2/', towriteID=False)


	with pytest.raises(TypeError):
		L = imgLoader(ra=ra , dec=dec, dir_obj=dir_obj, obsobj=obj, img_width=img_width, img_height=img_height)


def test_instantiate_imgLoader_error():
	"""
	test that an error being raised when none of obsobj or ra/dec/dir_obj are fed to imgLoader
	"""

	with pytest.raises(TypeError):
		L = imgLoader(img_width=img_width, img_height=img_height)


def test_call_download_stamp():
	"""
	test that download_stamp overwrites imgLoader
	"""

	newwidth = 5*u.arcsec
	newheight = 5*u.arcsec

	tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	obj = obsobj(tab, catalog='SDSS', dir_parent='./test2/', towriteID=False)

	L = imgLoader(obsobj=obj, img_width=img_width, img_height=img_height)

	L.download_stamp(band='r', img_width=newwidth, img_height=newheight)

	assert isinstance(L, imgLoader)
	assert L.ra == ra
	assert L.ra == obj.ra
	assert L.img_height == newheight
	assert L.img_width == newwidth
	assert L.band == 'r'
	assert L.dir_obj == obj.dir_obj
