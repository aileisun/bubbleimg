
import pytest
import shutil
import os
import numpy as np

from astropy.io import fits
import astropy.convolution as ac


from ...obsobj import obsObj

from .. import decomposer

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
def decomposer1(obj_dirobj):
	obj = obj_dirobj
	return decomposer.Decomposer(obj=obj, survey='hsc', z=z)


def test_decomposerr_init_getsurvey_fromobj(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj
	obj.z = z
	obj.survey = 'cfht'

	d = decomposer.Decomposer(obj=obj)
	assert d.survey == 'cfht'

	d = decomposer.Decomposer(obj=obj, survey='hsc')
	assert d.survey == 'hsc'


def test_spector_init_error_no_survey(obj_dirobj):
	""" complain when no survey is speficied either through obj or as attr of spector"""
	obj = obj_dirobj

	with pytest.raises(Exception) as e:
		d = decomposer.Decomposer(obj=obj)


def test_decomposer_get_z_fromobj(obj_dirobj):
	obj = obj_dirobj
	obj.add_sdss()

	d = decomposer.Decomposer(obj=obj, survey='hsc')
	assert d.z == z

	d = decomposer.Decomposer(obj=obj, survey='hsc', z=0.)
	assert d.z == 0.


def test_decomposer_get_conti_fnu_ratio_from_spector(decomposer1):
	d = decomposer1

	b1 = 'i'
	b2 = 'z'
	
	x = d._get_conti_fnu_ratio_from_spector(band1=b1, band2=b2)
	assert round(x, 3) == round(0.7302813376166697, 3)


