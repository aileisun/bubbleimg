import pytest
import shutil
import os
import astropy.table as at
import astropy.units as u
from collections import OrderedDict
from astropy.io import fits
import numpy as np

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

	fp = s.get_fp_stamp_noised(imgtag, img_sigma)

	assert os.path.isfile(fp)


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


def test_simulator_summarize_sim_noised(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	img_sigma=1
	niter = 5

	msrkwargs = OrderedDict(isocut=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'), minarea=5, onlycenter=True, centerradius=5.*u.arcsec)
	
	s.sim_noised(imgtag=imgtag, img_sigma=img_sigma, niter=niter, msrtype='iso', running_indx=False, keep_img=False, summarize=True, overwrite=True, **msrkwargs)

	noisedtag = s.get_tag_noised(img_sigma) # '_noised-1.0'
	fn_smr = s.dir_obj+'msr_iso{}_smr.csv'.format(noisedtag)
	assert os.path.isfile(fn_smr)
	tab = at.Table.read(fn_smr)
	assert 'area_ars_mean' in tab.colnames


def test_simulator_summarize_sim_noised_stand_alone(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	img_sigma=1
	niter = 5

	msrkwargs = OrderedDict(isocut=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'), minarea=5, onlycenter=True, centerradius=5.*u.arcsec)
	
	s.sim_noised(imgtag=imgtag, img_sigma=img_sigma, niter=niter, msrtype='iso', running_indx=False, keep_img=False, summarize=False, overwrite=True, **msrkwargs)

	s.summarize_sim_noised(imgtag=imgtag, img_sigma=img_sigma, msrtype='iso', overwrite=True, **msrkwargs)

	noisedtag = s.get_tag_noised(img_sigma) # '_noised-1.0'
	fn_smr = s.dir_obj+'msr_iso{}_smr.csv'.format(noisedtag)
	assert os.path.isfile(fn_smr)
	tab = at.Table.read(fn_smr)
	assert 'area_ars_mean' in tab.colnames	


def test_simulator_make_binned(simulator1):
	s = simulator1
	imgtag = 'OIII5008_I'
	binsize = 2
	s.make_binned(imgtag=imgtag, binsize=binsize, binpsf=True, overwrite=True)

	#=== image 
	fp = s.get_fp_stamp_binned(imgtag, binsize=binsize)
	assert os.path.isfile(fp)

	img = fits.getdata(s.get_fp_stamp(imgtag))
	img_b = fits.getdata(fp)
	assert np.all(np.array([img.shape])//binsize == np.array([img_b.shape]))
	assert np.sum(img) == np.sum(img_b) * 4. 

	#=== psf
	fp = s.get_fp_psf_binned(imgtag, binsize=binsize)
	assert os.path.isfile(fp)

	psf = fits.getdata(s.get_fp_psf(imgtag))
	psf_b = fits.getdata(fp)

	assert np.all(np.array([psf.shape])//binsize == np.array([psf_b.shape]))
	assert np.absolute(np.sum(psf_b) * 4. - np.sum(psf))/np.sum(psf) < 1.e-3


def test_simulator_make_smeared(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	gamma = 2.
	alpha = 1.

	s.make_smeared(imgtag=imgtag, gamma=gamma, alpha=alpha, overwrite=True)

	smearedtag = s.get_tag_smeared(gamma=gamma, alpha=alpha)

	fn = s.dir_obj+'stamp-OIII5008_I{}.fits'.format(smearedtag)
	assert os.path.isfile(fn)

	fn = s.dir_obj+'psf-OIII5008_I{}.fits'.format(smearedtag)
	assert os.path.isfile(fn)

	d = s._get_decomposer()
	d.make_psf_tab(imgtag=imgtag, msrsuffix='')
	d.make_psf_tab(imgtag=imgtag+smearedtag, msrsuffix='_smeared')

	fn = s.dir_obj+'psf.csv'
	tab = at.Table.read(fn)

	fn = s.dir_obj+'psf_smeared.csv'
	assert os.path.isfile(fn)
	tab_smeared = at.Table.read(fn)
	assert tab_smeared['psf_fwhm_pix'] > tab['psf_fwhm_pix']



def test_simulator_sim_smeared(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	smearargs = at.Table([[1.5, 2., 3., 4.], [2.5, 2.5, 2.5, 2.5]], names=['gamma', 'alpha'])
	s.sim_smeared(imgtag=imgtag, smearargs=smearargs, msrtype='iso', keep_img=True, overwrite=True)
	d = s._get_decomposer()
	d.make_psf_tab(imgtag=imgtag)

	for i, smeararg in enumerate(smearargs):
		alpha = smeararg['alpha']
		gamma = smeararg['gamma']

		smearedtag = s.get_tag_smeared(gamma=gamma, alpha=alpha)

		fn = s.dir_obj+'stamp-OIII5008_I{}.fits'.format(smearedtag)
		assert os.path.isfile(fn)

		fn = s.dir_obj+'psf-OIII5008_I{}.fits'.format(smearedtag)
		assert os.path.isfile(fn)

		fn = s.dir_obj+'msr_iso_smeared.csv'
		assert os.path.isfile(fn)
		tab = at.Table.read(fn)
		assert tab['imgtag'][i] == imgtag+smearedtag

		fn = s.dir_obj+'psf_smeared.csv'
		assert os.path.isfile(fn)


def test_simulator_sim_smeared_not_keep_img(simulator1):
	s = simulator1

	imgtag = 'OIII5008_I'
	smearargs = at.Table([[1.5, 2., 3., 4.], [2.5, 2.5, 2.5, 2.5]], names=['gamma', 'alpha'])
	s.sim_smeared(imgtag=imgtag, smearargs=smearargs, msrtype='iso', keep_img=False, overwrite=True)
	d = s._get_decomposer()
	d.make_psf_tab(imgtag=imgtag)

	for i, smeararg in enumerate(smearargs):
		alpha = smeararg['alpha']
		gamma = smeararg['gamma']

		smearedtag = s.get_tag_smeared(gamma=gamma, alpha=alpha)

		fn = s.dir_obj+'stamp-OIII5008_I{}.fits'.format(smearedtag)
		assert not os.path.isfile(fn)

		fn = s.dir_obj+'psf-OIII5008_I{}.fits'.format(smearedtag)
		assert not os.path.isfile(fn)

		fn = s.dir_obj+'msr_iso_smeared.csv'
		assert os.path.isfile(fn)
		tab = at.Table.read(fn)
		assert tab['imgtag'][i] == imgtag+smearedtag

		fn = s.dir_obj+'psf_smeared.csv'
		assert os.path.isfile(fn)
