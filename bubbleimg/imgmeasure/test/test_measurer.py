
import pytest
import shutil
import os
import astropy.units as u

from ...obsobj import obsObj

from .. import measurer

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'test_verification_data/SDSSJ0920+0034/'
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

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


@pytest.fixture
def measurer1(obj_dirobj):
	obj = obj_dirobj
	return measurer.Measurer(obj=obj, survey='hsc', z=z)


def test_measurer_init_getsurvey_fromobj(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj
	obj.z = z
	obj.survey = 'sdss'

	d = measurer.Measurer(obj=obj)
	assert d.survey == 'sdss'

	d = measurer.Measurer(obj=obj, survey='hsc')
	assert d.survey == 'hsc'


def test_measurer_init_error_no_survey(obj_dirobj):
	""" complain when no survey is speficied either through obj or as attr of measurer"""
	obj = obj_dirobj

	with pytest.raises(Exception) as e:
		d = measurer.Measurer(obj=obj)


def test_measurer_get_z_fromobj(obj_dirobj):
	obj = obj_dirobj
	obj.add_sdss()

	d = measurer.Measurer(obj=obj, survey='hsc')
	assert d.z == z

	d = measurer.Measurer(obj=obj, survey='hsc', z=0.)
	assert d.z == 0.


def test_measurer_get_fp_msr(measurer1):
	m = measurer1

	assert m.get_fp_stamp_img(imgtag='OIII5008_I') == m.dir_obj+'stamp-OIII5008_I.fits'

	assert m.get_fp_stamp_line_I(line='OIII5008') == m.dir_obj+'stamp-OIII5008_I.fits'


def test_measurer_theta_to_pix(measurer1):
	m = measurer1

	m._theta_to_pix(0.168*u.arcsec) == 1.
	m._theta_to_pix(2.*0.168*u.arcsec) == 2.


def test_measurer_pix_to_theta(measurer1):
	m = measurer1

	m._pix_to_theta(1.) == 0.168*u.arcsec
	m._pix_to_theta(1., wunit=False) == 0.168
