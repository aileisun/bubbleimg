import pytest
import shutil
import os
# import astropy.units as u

from ...obsobj import obsObj

from .. import simulator
from ...imgmeasure import isoMeasurer

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'verification_data/SDSSJ0920+0034/'
# bandline = 'i'
# bandconti = 'r'
# survey = 'hsc'
# z = 0.4114188

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
def simulator1(obj_dirobj):
	obj = obj_dirobj
	return simulator.Simulator(obj=obj, survey='hsc') # , z=z)


@pytest.fixture
def measurer1(obj_dirobj):
	obj = obj_dirobj
	return isoMeasurer(obj=obj, survey='hsc', z=z) # , z=z)


def test_simulator_init_getsurvey_fromobj(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj
	# obj.z = z
	obj.survey = 'sdss'

	s = simulator.Simulator(obj=obj)
	assert s.survey == 'sdss'

	s = simulator.Simulator(obj=obj, survey='hsc')
	assert s.survey == 'hsc'


def test_simulator_add_noise(simulator1):
	s = simulator1
	imgtag = 'OIII5008_I'
	img_sigma=1
	s.make_noised(imgtag=imgtag, img_sigma=img_sigma, overwrite=True)

	fp = s.get_fp_noised(imgtag, img_sigma)

	assert os.path.isfile(fp)

