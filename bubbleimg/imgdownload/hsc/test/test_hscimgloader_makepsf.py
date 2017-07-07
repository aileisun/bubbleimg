# test_hscimgloader.py
# ALS 2017/05/02

"""
to be used with pytest

test sets for hscimgloader

"""
import numpy as np
import astropy.table as at
import astropy.units as u
import shutil
import os
import pytest
from astropy.io import fits
import filecmp
import glob

from ..hscimgloader import hscimgLoader

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


def test_downlaod_psf(L_radec):
	""" test it can make_psf"""
	L = L_radec
	L._download_psf(band = 'r')

	assert os.path.isfile(L.dir_obj+'psf-r.fits')


def test_make_psf(L_radec):
	""" test it can make_psf"""
	L = L_radec
	L.make_psf(band = 'r', overwrite=True)

	assert os.path.isfile(L.dir_obj+'psf-r.fits')


def test_make_psf_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	assert L.make_psf(band = 'r', overwrite=True) is False


def test_make_psfs(L_radec):
	""" test it can make_psf"""
	L = L_radec
	L.make_psfs(overwrite=True)

	for b in ['g', 'r', 'i', 'z', 'y']:
		assert os.path.isfile(L.dir_obj+'psf-{0}.fits'.format(b))


def test_make_psfs_check_distinct(L_radec):
	""" test it can make_psf"""
	L = L_radec
	L.make_psfs(overwrite=True)

	for b in ['g', 'r', 'i', 'z', 'y']:
		assert os.path.isfile(L.dir_obj+'psf-{0}.fits'.format(b))

	for b0 in ['g', 'r', 'i', 'z', 'y']:
		for b1 in ['g', 'r', 'i', 'z', 'y']:
			if b0 != b1:
				psf0 = fits.getdata(L.get_fp_psf(b0))
				psf1 = fits.getdata(L.get_fp_psf(b1))

				assert not np.all(psf0 == psf1)



def test_make_psfs_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	assert L.make_psfs(overwrite=True) is False


def test_make_psfs_overwriteTrue(L_radec):
	"""
	test that when overwrite=True make_psf() always call download_psf() whether file exists or not
	"""
	L = L_radec

	band = 'r'
	overwrite = True

	file = dir_obj+'psf-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it creates psf
	assert not os.path.isfile(file)
	L.make_psfs(overwrite=overwrite)
	assert os.path.isfile(file)

	# when file does exist it overwrites psf
	if os.path.isfile(file):
		os.remove(file)
	open(file, 'w').close()

	L.make_psfs(overwrite=overwrite)
	assert os.path.isfile(file)
	assert os.stat(file).st_size > 0


def test_make_a_zero_size_file():
	file = 'testfile.txt'
	open(file, 'w').close()
	assert os.stat(file).st_size == 0
	os.remove(file)


def test_make_psfs_overwriteFalse(L_radec):

	"""
	test that when overwrite=False make_psf() does not update file
	"""
	L = L_radec

	overwrite = False

 	for band in L.bands:
		file = dir_obj+'psf-{0}.fits'.format(band)

		if not os.path.isdir(dir_obj):
			os.makedirs(dir_obj)

		if os.path.isfile(file):
			os.remove(file)

		open(file, 'w').close()
		assert os.stat(file).st_size == 0

	# when file exists it should not update file
	L.make_psfs(overwrite=overwrite)

 	for band in L.bands:
		file = dir_obj+'psf-{0}.fits'.format(band)
		assert os.stat(file).st_size == 0


# def test_make_psf_correctimgsize():

# 	for pixnum in [64, 128, 256]:
# 		L = hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=pixnum, img_height=pixnum)
# 		L.make_psf(band='r', overwrite=True)
# 		data = fits.getdata(L.dir_obj+'psf-r.fits')
# 		assert data.shape == (pixnum, pixnum)


# def test_make_psf_correctcontent():
# 	L = hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=128, img_height=128)

# 	L.make_psfs(overwrite=True)

#  	for band in L.bands:
#  		f = 'psf-{0}.fits'.format(band)
# 		file_totest = dir_obj+f
# 		file_verification = './test_verification_data_128pix/SDSSJ0920+0034/'+f
# 		assert filecmp.cmp(file_totest, file_verification)

# 	f = 'hsc_xid.csv'
# 	tab = at.Table.read(f, format='ascii.csv')
# 	for col in ['ra', 'dec', 'patch_id', 'tract', 'patch', 'patch_s', 'parent_id', ]:
# 		assert col in tab.colnames


