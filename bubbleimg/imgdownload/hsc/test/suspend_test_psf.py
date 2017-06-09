""" currently suspended as verfication data calexp-r is too big """

import pytest
from astropy.io import fits
import os 
from .. import psf



obj_id = 'SDSSJ0920+0034' 
ra     = 140.099364238908123 
dec    = 0.580160150759375104
tract  = 9564
patch_s  = '7,3'
fn_out = 'psf-g.fits'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./test/ and ./test2/ before and after test"""

	# setup
	if os.path.isfile(fn_out):
		os.remove(fn_out)
		
	yield
	# tear down
	if os.path.isfile(fn_out):
		os.remove(fn_out)


def test_psf():
	fn_in = 'test_verification_data_128pix/SDSSJ0920+0034/calexp-r.fits'
	psf.exposureF_to_psf(fn_in=fn_in, fn_out=fn_out, ra=ra, dec=dec)
	assert os.path.isfile(fn_out)