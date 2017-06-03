# test_loader_hsc.py
# ALS 2017/05/02

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


def test_make_stamp(L_radec):
	""" test it can make_stamp"""
	L = L_radec
	L.make_stamp(band = 'r', overwrite=True)
	L.make_stamp(band = 'i', overwrite=True)
	L.make_stamp(band = 'z', overwrite=True)

	L.plot_colorimg(img_type='stamp', bands ='riz')

	assert os.path.isfile(L.dir_obj+'color_stamp-riz.png')

