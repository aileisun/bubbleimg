import pytest
import shutil
import os
import astropy.table as at
import astropy.units as u

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

@pytest.fixture(scope="function", autouse=True)
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


def test_simulator_make_noised(simulator1):
	s = simulator1
	imgtag = 'OIII5008_I'
	img_sigma=1
	s.make_noised(imgtag=imgtag, img_sigma=img_sigma, overwrite=True)

	fp = s.get_fp_noised(imgtag, img_sigma)

	assert os.path.isfile(fp)


def test_simulator_get_measurer(simulator1):
	s = simulator1
	m = s.get_measurer()

	assert m.z == s.z


def test_simulator_sim_noised_keep_img(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	img_sigma=1
	niter = 5
	s.sim_noised(imgtag=imgtag, img_sigma=img_sigma, niter=niter, msrtype='iso', running_indx=True, keep_img=True, overwrite=True)

	noisedtag = s.get_tag_noised(img_sigma)
	fn = s.dir_obj+'msr_iso{}.csv'.format(noisedtag)

	assert os.path.isfile(s.get_fp_stamp_img(imgtag=imgtag+noisedtag+'_1'))

	assert os.path.isfile(fn)
	tab = at.Table.read(fn)

	assert len(tab) == niter


def test_simulator_sim_noised(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	img_sigma=1
	niter = 20
	s.sim_noised(imgtag=imgtag, img_sigma=img_sigma, niter=niter, msrtype='iso', running_indx=False, keep_img=False, overwrite=True)

	noisedtag = s.get_tag_noised(img_sigma)
	fn = s.dir_obj+'msr_iso{}.csv'.format(noisedtag)

	assert not os.path.isfile(s.get_fp_stamp_img(imgtag=imgtag+noisedtag+'_1'))

	assert os.path.isfile(fn)
	tab = at.Table.read(fn)

	assert len(tab) == niter



def test_simulator_sim_noised_customized_msr(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	img_sigma=1
	niter = 5

	msrkwargs = dict(isocut=1.e-15*u.Unit('erg / (arcsec2 cm2 s)'), minarea=20, onlycenter=True, centerradius=5.*u.arcsec)
	
	s.sim_noised(imgtag=imgtag, img_sigma=img_sigma, niter=niter, msrtype='iso', running_indx=False, keep_img=False, overwrite=True, **msrkwargs)

	noisedtag = s.get_tag_noised(img_sigma)
	fn = s.dir_obj+'msr_iso{}.csv'.format(noisedtag)

	assert not os.path.isfile(s.get_fp_stamp_img(imgtag=imgtag+noisedtag+'_1'))

	assert os.path.isfile(fn)
	tab = at.Table.read(fn)

	assert len(tab) == niter

	for key in msrkwargs:
		assert str(tab[key][0]) == str(msrkwargs[key])
