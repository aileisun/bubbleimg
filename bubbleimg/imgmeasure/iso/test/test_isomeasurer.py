
import pytest
import shutil
import os
import numpy as np
import astropy.units as u
import astropy.table as at
from astropy.io import fits
import filecmp

from ....obsobj import obsObj

from .. import isoMeasurer
from .. import polytools

dir_parent = './testing/'

bandline = 'i'
bandconti = 'r'
survey = 'hsc'

isocut = 1.e-15*u.Unit('erg / (arcsec2 cm2 s)')

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	os.makedirs(dir_parent)

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)



@pytest.fixture
def measurer1():
	ra = 140.099341430207
	dec = 0.580162492432517
	z = 0.4114188

	dir_obj = './testing/SDSSJ0920+0034/'
	dir_verif = 'verification_data/SDSSJ0920+0034/'

	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)

	shutil.copytree(dir_verif, dir_obj)
	obj = obsObj(ra=ra, dec=dec, dir_obj = dir_obj)
	return isoMeasurer(obj=obj, survey='hsc', z=z, center_mode='n/2-1')


@pytest.fixture
def measurer_nanimg():
	ra = 140.826752453534
	dec = 0.530545728517824

	z = 0.548
	dir_verif = './verification_data/SDSSJ092318+003149/'
	dir_obj = './testing/SDSSJ092318+003149/'

	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)
	shutil.copytree(dir_verif, dir_obj)

	obj = obsObj(ra=ra, dec=dec, dir_obj=dir_obj, obj_naming_sys='sdss_precise')
	return isoMeasurer(obj=obj, survey='hsc', z=z, center_mode='n/2-1')


def test_isomeasurer_get_fp_msr(measurer1):
	m = measurer1

	assert m.get_fp_msr(imgtag='OIII5008_I') == m.dir_obj+'msr_iso-OIII5008_I.csv'


def test_isomeasurer_get_fp_contours(measurer1):
	m = measurer1


	imgtag = 'OIII5008_I'

	assert m.get_fp_contours(imgtag=imgtag, onlycenter=False) == m.dir_obj+'msr_iso-OIII5008_I_contours.pkl'
	assert m.get_fp_contours(imgtag=imgtag, onlycenter=True) == m.dir_obj+'msr_iso-OIII5008_I_contours-ctr.pkl'


def test_isomeasurer_get_contours_from_img(measurer1):
	m = measurer1

	fn_img = m.dir_obj+'stamp-OIII5008_I.fits'
	img = fits.getdata(fn_img)

	xc, yc = m._get_xc_yc(img)
	isocut = 1.e-15
	contours = m._get_contours_from_img(img=img, isocut=isocut, xc=xc, yc=yc, minarea=0., onlycenter=False, centerradius=2.*u.arcsec)

	assert polytools.ispolys(contours)
	assert polytools.NetPolygonsArea(contours) > 0


def test_isomeasurer_make_measurements(measurer1):
	m = measurer1
 	imgtag = 'OIII5008_I'
	minarea = 5
	status = m.make_measurements(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=True, centerradius=5.*u.arcsec, overwrite=False, savecontours=True, plotmsr=True)

	assert status

	assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I.pdf')
	assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I.csv')
	assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I_contours-ctr.pkl')

	tab = at.Table.read(m.dir_obj+'msr_iso-OIII5008_I.csv')

	assert tab['imgtag'][0] == imgtag
	assert u.Quantity(tab['isocut'][0]) == isocut
	assert tab['minarea'][0] == minarea
	assert tab['area_ars'] == tab['area_pix'] * 0.168**2
	assert tab['dmax_ars'] == tab['dmax_pix'] * 0.168

	assert round(tab['area_kpc'], 1) == round(tab['area_ars'] * 5.465**2, 1)
	assert round(tab['dmax_kpc'], 1) == round(tab['dmax_ars'] * 5.465, 1)



def test_isomeasurer_make_measurements_suffix(measurer1):
	m = measurer1
 	imgtag = 'OIII5008_I'
	minarea = 5
	isocut1 = 1.e-15*u.Unit('erg / (arcsec2 cm2 s)')
	suffix1 = '_1e-15'

	isocut2 = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')
	suffix2 = '_3e-15'

	for isocut, suffix in ((isocut1, suffix1), (isocut2, suffix2)):
		status = m.make_measurements(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=True, centerradius=5.*u.arcsec, suffix=suffix, overwrite=True, savecontours=True, plotmsr=True)

		assert status

		assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I{suffix}.pdf'.format(suffix=suffix))
		assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I{suffix}.csv'.format(suffix=suffix))
		assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I{suffix}_contours-ctr.pkl'.format(suffix=suffix))


	f1 = m.dir_obj+'msr_iso-OIII5008_I{suffix}.csv'.format(suffix=suffix1)
	f2 = m.dir_obj+'msr_iso-OIII5008_I{suffix}.csv'.format(suffix=suffix2)
	assert not filecmp.cmp(f1, f2)



def test_isomeasurer_nan_image(measurer_nanimg):
	m = measurer_nanimg

 	imgtag = 'OIII5008_I'
	minarea = 5
	status = m.make_measurements(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=True, centerradius=5.*u.arcsec, overwrite=False, savecontours=True, plotmsr=True)

	fn = m.get_fp_msr(imgtag=imgtag)
	tab = at.Table.read(fn)

	for col in ['area_kpc', 'dmax_kpc', 'rmax_kpc', 'dper_kpc', 'area_ars', 'dmax_ars', 'rmax_ars', 'dper_ars', 'area_pix', 'dmax_pix', 'rmax_pix', 'dper_pix', 'theta_dmax', 'theta_rmax', 'theta_dper', 'aspectr', ]:
		assert np.isnan(tab[col][0])

	for col in ['imgtag', 'isocut', 'minarea', 'onlycenter', 'centerradius']:
		assert col in tab.colnames


def test_isomeasurer_make_colorimg(measurer1):
	m = measurer1

	status = m.make_colorimg(bands ='gri', img_type='stamp', overwrite=False)

	assert status
	assert os.path.isfile(m.dir_obj+'color_stamp-gri.png')


def test_isomeasurer_make_visualpanel(measurer1):
	m = measurer1

	# isocut = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')

	status = m.make_visualpanel(compo_bands ='riz', imgtag='OIII5008_I', overwrite=False)

	assert status
	assert os.path.isfile(m.dir_obj+'msr_iso-OIII5008_I_panel.pdf')


def test_isomeasurer_make_noiselevel(measurer1):
	m = measurer1

	# isocut = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')

	status = m.make_noiselevel(imgtag='OIII5008_I', toplot=True, overwrite=False)

	assert status
	assert os.path.isfile(m.dir_obj+'noiselevel-OIII5008_I.csv')
	assert os.path.isfile(m.dir_obj+'noiselevel-OIII5008_I.pdf')

	assert m.get_noiselevel(imgtag='OIII5008_I', wunit=False) > 0.
	assert m.get_noiselevel(imgtag='OIII5008_I', wunit=True).unit == u.Unit('1e-15 erg / (arcsec2 cm2 s)')


