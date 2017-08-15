
import pytest
import shutil
import os
import numpy as np
import astropy.table as at


from astropy.io import fits
import astropy.convolution as ac
import astropy.units as u

from ....obsobj import obsObj

from .. import plainDecomposer

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'verification_data/SDSSJ0920+0034/'
bandline = 'i'
bandconti = 'r'
survey = 'hsc'
z = 0.4114188

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	os.makedirs(dir_parent)
	shutil.copytree(dir_verif, dir_obj)

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


@pytest.fixture
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


@pytest.fixture
def decomposer1(obj_dirobj):
	obj = obj_dirobj
	return plainDecomposer(obj=obj, survey='hsc', z=z)


def test_plaindecomposer_get_fp_psf(decomposer1):
	d = decomposer1

	fp_stamp = d.get_fp_stamp(band='i')
	fp_psf = d.get_fp_psf(fp_stamp)
	print fp_stamp
	print fp_psf

	assert os.path.dirname(fp_stamp) == os.path.dirname(fp_psf)
	assert os.path.basename(fp_stamp)[5:] == os.path.basename(fp_psf)[3:]
	assert os.path.basename(fp_psf)[:3] == 'psf'

	assert 'psf-oiii5008_over-halpha' == d.get_fp_psf('stamp-oiii5008_over-halpha')


def test_decomposer_init_getsurvey_fromobj(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj
	obj.z = z
	obj.survey = 'sdss'

	d = plainDecomposer(obj=obj)
	assert d.survey == 'sdss'

	d = plainDecomposer(obj=obj, survey='hsc')
	assert d.survey == 'hsc'


def test_spector_init_error_no_survey(obj_dirobj):
	""" complain when no survey is speficied either through obj or as attr of spector"""
	obj = obj_dirobj

	with pytest.raises(Exception) as e:
		d = plainDecomposer(obj=obj)


def test_decomposer_get_z_fromobj(obj_dirobj):
	obj = obj_dirobj
	obj.add_sdss()

	d = plainDecomposer(obj=obj, survey='hsc')
	assert d.z == z

	d = plainDecomposer(obj=obj, survey='hsc', z=0.)
	assert d.z == 0.


def test_decomposer_get_conti_fnu_ratio_from_spector(decomposer1):
	d = decomposer1

	b1 = 'i'
	b2 = 'z'
	
	x = d._get_conti_fnu_ratio_from_spector(band1=b1, band2=b2)
	assert round(x, 3) == round(0.7302813376166697, 3)



def test_decomposer_make_stamp_psfmatch(decomposer1):
	d = decomposer1

	band = bandline
	bandto = bandconti


	status = d.make_stamp_psfmatch(band, bandto, overwrite=True)

	assert status

	status = d.make_stamp_psfmatch(band, bandto, overwrite=False, )

	assert status

	fp = d.get_fp_stamp_psfmatched(band, bandto)

	for fn in [fp, d.get_fp_psf(fp)]:
		assert os.path.isfile(fn)

	status = d.make_stamp_psfmatch(band, bandto, overwrite=True, towrite_psk=True)

	for fn in [fp, d.get_fp_psf(fp), d.get_fp_psk(fp)]:
		assert os.path.isfile(fn)


def test_decomposer_subtract_img_w_ratio(decomposer1):
	d = decomposer1

	fn1 = d.dir_obj+'a1.fits'
	fn2 = d.dir_obj+'a2.fits'
	fnout = d.dir_obj+'answer.fits'

	arr1 = np.ones([2, 2])
	arr2 = np.array([[1, 0], [1, 0]])
	arr_ans = np.array([[1, 2], [1, 2]])

	fits.PrimaryHDU(arr1).writeto(fn1, overwrite=True)
	fits.PrimaryHDU(arr2).writeto(fn2, overwrite=True)

	d._subtract_img_w_ratio(fn1, fn2, fnout, a1=2., a2=1., overwrite=True)

	arr_out = fits.getdata(fnout)

	assert np.all(arr_out == arr_ans)


def test_decomposer_has_smaller_psf(decomposer1):
	d = decomposer1


	assert d._band_has_smaller_psf(band='i', bandto='r')
	assert d._band_has_smaller_psf(band='i', bandto='y')



def test_decomposer_make_stamp_contsub(decomposer1):

	d = decomposer1

	for bandconti in ['r', 'y']:
		fn = d.get_fp_stamp_contsub(bandline, bandconti)
		fn_psf = d.get_fp_psf(fn)

		status = d.make_stamp_contsub(bandline, bandconti, overwrite=True)

		assert status

		assert os.path.isfile(fn)
		assert os.path.isfile(fn_psf)


def test_decomposer_make_linemap(decomposer1):
	d = decomposer1

	bandline = 'i'

	for bandconti in ['r', 'y']:
		d.make_stamp_linemap(bandline=bandline, bandconti=bandconti, line='OIII5008', overwrite=True)


		assert os.path.isfile(d.dir_obj+'stamp-OIII5008.fits')
		assert os.path.isfile(d.dir_obj+'psf-OIII5008.fits')

		hdus = fits.open(d.dir_obj+'stamp-OIII5008.fits')
		assert u.Unit(hdus[0].header['BUNIT']) == u.Unit('1e-17 erg s-1 cm-2')

		fn_contsub = d.get_fp_stamp_contsub(bandline, bandconti)
		fn_psf_contsub = d.get_fp_psf(fn_contsub)
		assert np.all(fits.getdata(d.dir_obj+'psf-OIII5008.fits') == fits.getdata(fn_psf_contsub))


def test_decomposer_make_linemap_I(decomposer1):
	d = decomposer1
	for bandconti in ['r', 'y']:
		d.make_stamp_linemap(bandline=bandline, bandconti=bandconti, line='OIII5008', overwrite=True)
		d.make_stamp_linemap_I(bandline=bandline, bandconti=bandconti, line='OIII5008', overwrite=True)

		assert os.path.isfile(d.dir_obj+'stamp-OIII5008_I.fits')
		assert os.path.isfile(d.dir_obj+'psf-OIII5008_I.fits')



def test_decomposer_plot_psfmatch(decomposer1):
	d = decomposer1

	band = bandline
	bandto = bandconti

	fn = d.dir_obj+'psf_testing.pdf'

	status = d.make_stamp_psfmatch(band, bandto, overwrite=True)
	assert status

	status = d.plot_psfmatch(band, bandto, fn=fn, matching=True, overwrite=True)
	assert status

	assert os.path.isfile(fn)

	status = d.plot_psfmatch(band='y', bandto='z', fn=fn, matching=True, overwrite=True)
	assert status == False



def test_decomposer_make_stamp_contsub_writetab(decomposer1):

	d = decomposer1

	for bandconti in ['r', 'y']:

		status = d.make_stamp_contsub(bandline, bandconti, overwrite=True, towrite_tab=True, toplot_psf=True)

		fn = d.get_fp_stamp_contsub(bandline, bandconti)
		fn_tab = d.get_fp_contsubtab(bandline, bandconti)

		assert status
		assert os.path.isfile(fn)
		assert os.path.isfile(fn_tab)

		tab = at.Table.read(fn_tab, format='ascii.csv')
		assert tab['band_psfm_from'][0] == bandline
		assert tab['band_psfm_to'][0] == bandconti
		assert tab['psf_fwhm'][0] < tab['psf_fwhm_conti'][0]

		assert os.path.isfile(d.dir_obj + 'contsub_psf-{}-{}.pdf'.format(bandline, bandconti))


def test_decomposer_make_psf_tab(decomposer1):

	d = decomposer1
	status = d.make_psf_tab(overwrite=True)

	assert status
	assert os.path.isfile(d.fp_psf_tab)

	tab = at.Table.read(d.fp_psf_tab, format='ascii.csv')

	assert tab['psf_fwhm_i'][0] < tab['psf_fwhm_y'][0]

	for band in d.bands:
		assert 'psf_fwhm_'+band in tab.colnames