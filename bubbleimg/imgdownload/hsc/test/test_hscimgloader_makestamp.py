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

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


@pytest.fixture
def L_radec():
	""" returns a imgLoader object initiated with the ra dec above"""
	return hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=img_width, img_height=img_height)



def test_downlaod_stamp(L_radec):
	""" test it can make_stamp"""
	L = L_radec
	L._download_stamp(band = 'r')

	assert os.path.isfile(L.dir_obj+'stamp-r.fits')


def test_make_stamp(L_radec):
	""" test it can make_stamp"""
	L = L_radec
	L.make_stamp(band = 'r', overwrite=True)

	assert os.path.isfile(L.dir_obj+'stamp-r.fits')


def test_make_stamp_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	assert L.make_stamp(band = 'r', overwrite=True) is False


def test_make_stamps(L_radec):
	""" test it can make_stamp"""
	L = L_radec
	L.make_stamps(overwrite=True)

	for b in ['g', 'r', 'i', 'z', 'y']:
		assert os.path.isfile(L.dir_obj+'stamp-{0}.fits'.format(b))


def test_make_stamps_no_img_in_hsc():
	""" try loading a img thats not in hsc """
	ra = 150.0547735
	dec = 12.7073027
	dir_obj = './testing/SDSSJ1000+1242/'

	L = hscimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, img_width=img_width, img_height=img_height)

	assert L.make_stamps(overwrite=True) is False


def test_make_stamps_overwriteTrue(L_radec):
	"""
	test that when overwrite=True make_stamp() always call download_stamp() whether file exists or not
	"""
	L = L_radec

	band = 'r'
	overwrite = True

	file = dir_obj+'stamp-{0}.fits'.format(band)

	if os.path.isfile(file):
		os.remove(file)

	# when file does not exist it creates stamp
	assert not os.path.isfile(file)
	L.make_stamps(overwrite=overwrite)
	assert os.path.isfile(file)

	# when file does exist it overwrites stamp
	if os.path.isfile(file):
		os.remove(file)
	open(file, 'w').close()

	L.make_stamps(overwrite=overwrite)
	assert os.path.isfile(file)
	assert os.stat(file).st_size > 0


def test_make_a_zero_size_file():
	file = 'testfile.txt'
	open(file, 'w').close()
	assert os.stat(file).st_size == 0
	os.remove(file)


def test_make_stamps_overwriteFalse(L_radec):

	"""
	test that when overwrite=False make_stamp() does not update file
	"""
	L = L_radec

	overwrite = False

 	for band in L.bands:
		file = dir_obj+'stamp-{0}.fits'.format(band)

		if not os.path.isdir(dir_obj):
			os.makedirs(dir_obj)

		if os.path.isfile(file):
			os.remove(file)

		open(file, 'w').close()
		assert os.stat(file).st_size == 0

	# when file exists it should not update file
	L.make_stamps(overwrite=overwrite)

 	for band in L.bands:
		file = dir_obj+'stamp-{0}.fits'.format(band)
		assert os.stat(file).st_size == 0


def test_make_stamp_correctimgsize():

	for pixnum in [64, 128, 256]:
		L = hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=pixnum, img_height=pixnum)
		L.make_stamp(band='r', overwrite=True)
		data = fits.getdata(L.dir_obj+'stamp-r.fits')
		assert data.shape == (pixnum, pixnum)


# def test_make_stamp_correctcontent():
# 	L = hscimgLoader(ra=ra , dec=dec, dir_parent=dir_parent, img_width=128, img_height=128)

# 	L.make_stamps(overwrite=True)

#  	for band in L.bands:
#  		f = 'stamp-{0}.fits'.format(band)
# 		file_totest = dir_obj+f
# 		file_verification = './test_verification_data_128pix/SDSSJ0920+0034/'+f
# 		assert filecmp.cmp(file_totest, file_verification)

# 	f = 'hsc_xid.csv'
# 	tab = at.Table.read(L.dir_obj+f, format='ascii.csv')
# 	for col in ['ra', 'dec', 'patch_id', 'tract', 'patch', 'patch_s', 'parent_id', ]:
# 		assert col in tab.colnames


def test_header_w_BUNIT(L_radec):
	L = L_radec
	L.make_stamp(band='r')

	header = fits.getheader(L.dir_obj+'stamp-r.fits')
	assert header['BUNIT'] == '1.58479740e-02 nanomaggy'
	# assert 'FLUXMAG0' not in header


def test_unit_conversion(L_radec):
	L = L_radec

	filein = dir_obj+'test_filein.fits'
	fileout = dir_obj+'test_fileout.fits'

	data = np.array([[63095734448.0194]])
	hdu = fits.PrimaryHDU(data)
	hdu.header['FLUXMAG0'] = 63095734448.0194
	hdu1 = fits.ImageHDU(data)
	hdus = fits.HDUList([hdu, hdu1])
	hdus.writeto(filein, overwrite=True)

	assert hdu.header['FLUXMAG0'] == 63095734448.0194

	L._write_fits_unit_converted_to_nanomaggy(filein, fileout)

	dataout = fits.getdata(fileout)
	headerout = fits.getheader(fileout)

	assert round(dataout[0, 0],3) == 10.**9
	assert headerout['BUNIT'] == 'nanomaggy'


def test_make_stamp_not_tokeepraw(L_radec):
	L = L_radec
	
	if os.path.isdir(L.dir_obj):
		shutil.rmtree(L.dir_obj)

	L.make_stamp(band='r', tokeepraw=False, overwrite=True)

	a = glob.glob(L.dir_obj+'cutout*.fits')
	assert len(a) == 0


def test_make_stamp_tokeepraw(L_radec):
	L = L_radec
	
	if os.path.isdir(L.dir_obj):
		shutil.rmtree(L.dir_obj)

	L.make_stamp(band='r', tokeepraw=True, overwrite=True)

	a = glob.glob(L.dir_obj+'cutout*.fits')
	assert len(a) > 0

