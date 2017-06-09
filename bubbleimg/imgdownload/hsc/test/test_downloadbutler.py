
import pytest
from astropy.io import fits
import os 

from .. import downloadbutler



obj_id = 'SDSSJ0920+0034' 
ra     = 140.099364238908123 
dec    = 0.580160150759375104
tract  = 9564
patch_s  = '7,3'
fileout = 'calexp-g.fits'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./test/ and ./test2/ before and after test"""

	# setup
	if os.path.isfile(fileout):
		os.remove(fileout)

	yield
	# tear down
	if os.path.isfile(fileout):
		os.remove(fileout)


def test_downloadbutler_download_iaa():
	downloadbutler_download(environment='iaa')


def test_downloadbutler_download_online():
	downloadbutler_download(environment='online')


def downloadbutler_download(environment='iaa'):
	b = downloadbutler.multiButler(environment=environment, release_version='dr1', semester='s16a', rerun='s16a_wide').butler

	with b:
		dataId = dict(tract=tract, patch_s=patch_s, filter="HSC-G")


		status = b.download_file(fileout, **dataId)

		assert status

		assert os.path.isfile(fileout)
		os.remove(fileout)
		assert not os.path.isfile(fileout)


def test_downloadbutler_download_psf_iaa():
	fileout = 'test_psf.fits'
	b = downloadbutler.multiButler(environment='iaa', release_version='dr1', semester='s16a', rerun='s16a_wide').butler

	with b:
		dataId = dict(tract=tract, patch_s=patch_s, filter="HSC-G")

		status = b.download_psf(fileout, ra=ra, dec=dec, **dataId)

		assert status

		assert os.path.isfile(fileout)
		os.remove(fileout)
		assert not os.path.isfile(fileout)
