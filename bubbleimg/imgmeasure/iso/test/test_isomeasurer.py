
import pytest
import shutil
import os
import numpy as np
import astropy.units as u
import astropy.table as at

from ....obsobj import obsObj

from .. import isoMeasurer
from .. import polytools

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'verification_data/SDSSJ0920+0034/'
bandline = 'i'
bandconti = 'r'
survey = 'hsc'
z = 0.4114188

isocut = 1.e-15*u.Unit('erg / (arcsec2 cm2 s)')

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
def measurer1(obj_dirobj):
	obj = obj_dirobj
	return isoMeasurer(obj=obj, survey='hsc', z=z, center_mode='n/2-1')


def test_isomeasurer_get_fp_msr(measurer1):
	m = measurer1

	assert m.get_fp_msr(imgtag='OIII5008_I') == m.dir_obj+'msr_iso-OIII5008_I.csv'


def test_isomeasurer_get_fp_contour(measurer1):
	m = measurer1


	imgtag='OIII5008_I'

	assert m.get_fp_contour(imgtag=imgtag, onlycenter=False) == m.dir_obj+'contour-OIII5008_I.pkl'
	assert m.get_fp_contour(imgtag=imgtag, onlycenter=True) == m.dir_obj+'contour_ctr-OIII5008_I.pkl'



# def test_isomeasurer_measure_line_I(measurer1):
# 	m = measurer1

# 	m.make_measurements_line_I(line='OIII5008', overwrite=True, isocut=isocut, isoareallimit=0, onlycenter=True)

# 	assert os.path.file(m.dir_obj+'msr_iso-OIII5008_I.csv')



def test_isomeasurer_get_contours_from_fits(measurer1):
	m = measurer1

	fn_img = m.dir_obj+'stamp-OIII5008_I.fits'
	contours = m._get_contours_from_fits(fn_img=fn_img, isocut=isocut, minarea=0., onlycenter=False, centerradius=2.*u.arcsec)

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
	assert os.path.isfile(m.dir_obj+'contours_ctr-OIII5008_I.pkl')

	tab = at.Table.read(m.dir_obj+'msr_iso-OIII5008_I.csv')

	assert tab['imgtag'][0] == imgtag
	assert u.Quantity(tab['isocut'][0]) == isocut
	assert tab['minarea'][0] == minarea
	assert tab['area_ars'] == tab['area_pix'] * 0.168**2
	assert tab['dmax_ars'] == tab['dmax_pix'] * 0.168

	assert round(tab['area_kpc'], 1) == round(tab['area_ars'] * 5.465**2, 1)
	assert round(tab['dmax_kpc'], 1) == round(tab['dmax_ars'] * 5.465, 1)





# def test_find_contours_one_w_unit():

# 	unit = u.Unit('erg s-1 cm-2')
# 	img = img_dot * unit
# 	threshold = thrshld * unit

# 	contours = polytools.find_contours(img, threshold, tocomplete=True)

# 	assert len(contours) == 1.

# 	print("real area: ", np.sum(img), "cntr area: ", polytools.NetPolygonsArea(contours))

# 	assert about_the_same_area(contours, np.sum(img), uncertainty=0.)
