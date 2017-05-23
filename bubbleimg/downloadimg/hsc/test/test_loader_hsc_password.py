# test_loader_hsc_password.py
"""
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

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './test/'
dir_obj = './test/SDSSJ0920+0034/'
img_width = 128*u.pix
img_height = 128*u.pix



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


def test_HSCimgLoader_user_password():
	""" test it can make_stamp"""
	
	L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height, user=user, password=password)
	L.make_stamp(band = 'r', overwrite=True)

	assert os.path.isfile(L.dir_obj+'stamp-r.fits')
