
import pytest
import shutil
import os
import astropy.units as u
import numpy as np
import astropy.table as at
import filecmp

from ...obsobj import obsObj

from .. import spector
from .. import inttools

ra = 140.099341430207
dec = 0.580162492432517
z = 0.4114
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'
dir_verif = 'test_verification_data/SDSSJ0920+0034/'
survey = 'hsc'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	os.makedirs(dir_obj)
	# shutil.copytree(dir_verif, dir_obj)
	shutil.copyfile(dir_verif+'spec.fits', dir_obj+'spec.fits')

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


@pytest.fixture
def spector1(obj_dirobj):
	obj = obj_dirobj
	return spector.Spector(obj=obj, survey_spec='boss', survey='hsc')


def test_spector_init(spector1):
	s = spector1

	spec, lcoord = s.get_spec_lcoord()

	assert len(spec) == len(lcoord)


def test_spector_init_error_survey_spec(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj

	with pytest.raises(Exception) as e: 
		s = spector.Spector(obj=obj, survey_spec='sdss')


def test_spector_init_autochoose_survey_spec(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj

	s = spector.Spector(obj=obj, survey='hsc')

	assert s.survey_spec == 'boss'


def test_spector_init_getsurvey_fromobj(obj_dirobj):
	""" complain when the survey_spec is set to a wrong value"""
	obj = obj_dirobj
	obj.survey = 'cfht'

	s = spector.Spector(obj=obj)

	assert s.survey == 'cfht'

	s = spector.Spector(obj=obj, survey='hsc')

	assert s.survey == 'hsc'


def test_spector_init_error_no_survey(obj_dirobj):
	""" complain when no survey is speficied either through obj or as attr of spector"""
	obj = obj_dirobj

	with pytest.raises(Exception) as e:
		s = spector.Spector(obj=obj)



def test_spector_make_spec_decomposed_ecsv(spector1):
	s = spector1
	s.make_spec_decomposed_ecsv(overwrite=True)

	assert os.path.isfile(s.dir_obj+'spec_decomposed.ecsv')

	tab = at.Table.read(s.dir_obj+'spec_decomposed.ecsv', format='ascii.ecsv')

	units = [tab[col].unit for col in ['spec', 'speccont', 'specline']]

	# assert all units are identical
	assert len(set(units)) == 1


def test_spector_make_spec_contextrp_ecsv(spector1):
	s = spector1
	fn = s.dir_obj+'spec_contextrp.ecsv'
	s.make_spec_contextrp_ecsv(overwrite=True)

	assert os.path.isfile(fn)

	tab = at.Table.read(fn, format='ascii.ecsv')

	assert tab['lcoord'].unit == u.AA
	assert tab['speccontextrp'].unit == u.Unit("erg / (Angstrom cm2 s)")


def test_spector_get_contline(spector1):
	s = spector1

	spec, lcoord = s.get_spec_lcoord_from_spectab(component='all')
	speccont, __ = s.get_spec_lcoord_from_spectab(component='cont')
	specline, __ = s.get_spec_lcoord_from_spectab(component='line')

	assert len(speccont) == len(lcoord)
	assert len(specline) == len(lcoord)

	assert max(speccont) <= max(spec)
	assert max(specline) <= max(spec)


def test_spector_plot_spec(spector1):
	s = spector1
	
	s.plot_spec(wfilters=True, wconti=True, wcontextrp=True, wline=True, overwrite=True)
	assert os.path.isfile(s.dir_obj+'spec.pdf')


def test_inttools_calc_Fnu_in_band(spector1):
	s = spector1
	spec, lcoord = s.get_spec_lcoord()
	trans, lcoord_trans = s.get_filter_trans(band='i')

	Fnu = inttools.calc_Fnu_in_band_from_fl(fl=spec, ls=lcoord, trans=trans, ls_trans=lcoord_trans)


	assert (Fnu/u.Unit("erg s-1 cm-2 Hz-1")).unit == u.dimensionless_unscaled



def test_spector_calc_Fnu_in_band(obj_dirobj):
	obj = obj_dirobj
	s = spector.Spector(obj=obj, survey_spec='boss', survey='sdss')
	s.obj.sdss.load_photoobj()

	for band in s.bands:
		Fnu = s.calc_Fnu_in_band(band=band, component='all')
		mAB = s.calc_mAB_in_band(band=band, component='all')

		convovledflux = Fnu.to(u.nanomaggy).value
		fiber2flux = s.obj.sdss.photoobj['fiber2Flux_'+band] # in nanomaggy
		fiber2mag = s.obj.sdss.photoobj['fiber2Mag_'+band]
		assert np.absolute(convovledflux - fiber2flux) < 10.
		assert np.absolute(mAB.value - fiber2mag) < 1.


def test_spector_make_spec_mag(spector1): 
	s = spector1
	s.make_spec_mag(overwrite=True)

	fn = s.dir_obj+'spec_mag.csv'
	assert os.path.isfile(fn)

	tab = at.Table.read(fn, comment='#', format='ascii.csv')

	for band in s.bands:
		assert 'specMag_{0}'.format(band) in tab.colnames
		assert 'speccontMag_{0}'.format(band) in tab.colnames
		assert 'speclineMag_{0}'.format(band) in tab.colnames
		assert 'speccontextrpMag_{0}'.format(band) in tab.colnames

		assert 'specFnu_{0}'.format(band) in tab.colnames
		assert 'speccontFnu_{0}'.format(band) in tab.colnames
		assert 'speclineFnu_{0}'.format(band) in tab.colnames
		assert 'speccontextrpFnu_{0}'.format(band) in tab.colnames



def test_spector_get_spec_mag(spector1): 
	s = spector1
	tab = s.get_spec_mag_tab()

	for band in s.bands:
		assert 'specMag_{0}'.format(band) in tab.colnames
		assert 'speccontMag_{0}'.format(band) in tab.colnames
		assert 'speclineMag_{0}'.format(band) in tab.colnames

		assert 'specFnu_{0}'.format(band) in tab.colnames
		assert 'speccontFnu_{0}'.format(band) in tab.colnames
		assert 'speclineFnu_{0}'.format(band) in tab.colnames



def test_spector_get_spec_mag_value(spector1):
	s = spector1

	x = s.get_spec_mag_value(component='line', fluxquantity='mag', band='i')
	assert isinstance(x, np.float64)
	assert round(x, 3) == round(21.7336867073556, 3)

	x = s.get_spec_mag_value(component='cont', fluxquantity='fnu', band='g')
	assert round(x, 3) == round(1.678181807435394, 3)


def test_spector_get_fnu_ratio_band1_over_band2(spector1):
	s = spector1
	b1 = 'i'
	b2 = 'z'
	x = s.get_fnu_ratio_band1_over_band2(band1=b1, band2=b2, component='cont')
	
	assert round(x, 3) == round(0.7302813376166697, 3)