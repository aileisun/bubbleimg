import shutil
import os
import pytest
import filecmp

from .. import make_psf

ra     = 140.099364238908123 
dec    = 0.580160150759375104
fn_in = "/array2/SSP/dr1/s16a/data/s16a_wide/deepCoadd/HSC-R/9564/7,3/calexp-HSC-R-9564-7,3.fits"

dir_testing = '/data/home/hscpipe/alsun/get_psf/test/testing/'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after test"""

	# setup
	if os.path.isdir(dir_testing):
		shutil.rmtree(dir_testing)

	os.makedirs(dir_testing)

	yield
	# tear down
	if os.path.isdir(dir_testing):
		shutil.rmtree(dir_testing)


def test_call_make_psf():

	fn_out = dir_testing+'test.fits'

	status = make_psf.main(fn_in, fn_out, ra, dec)

	assert status
	assert os.path.isfile(fn_out)
	if os.path.isfile(fn_out): 
		os.remove(fn_out)


def test_call_make_distinct_psf():

	fn_in_r = "/array2/SSP/dr1/s16a/data/s16a_wide/deepCoadd/HSC-R/9564/7,3/calexp-HSC-R-9564-7,3.fits"
	fn_out_r = dir_testing+'test_r.fits'
	status_r = make_psf.main(fn_in_r, fn_out_r, ra, dec)

	assert status_r
	assert os.path.isfile(fn_out_r)


	fn_in_z = "/array2/SSP/dr1/s16a/data/s16a_wide/deepCoadd/HSC-Z/9564/7,3/calexp-HSC-Z-9564-7,3.fits"
	fn_out_z = dir_testing+'test_z.fits'
	status_z = make_psf.main(fn_in_z, fn_out_z, ra, dec)

	assert status_z
	assert os.path.isfile(fn_out_z)

	# the two files are distinct
	assert not filecmp.cmp(fn_out_r, fn_out_z)

