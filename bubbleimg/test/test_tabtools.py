
import pytest
import os
import shutil
import astropy.table as at
from astropy.io import ascii
import bubbleimg
from astropy import units as u
import numpy as np

from .. import tabtools


dir_test = './testing/'
dir_verif = 'verification_data_tabtools/'

@pytest.fixture(scope="function", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_test):
		shutil.rmtree(dir_test)

	shutil.copytree(dir_verif, dir_test)

	yield
	# tear down
	if os.path.isdir(dir_test):
		shutil.rmtree(dir_test)


@pytest.fixture
def measurer1():
	ra = 140.099341430207
	dec = 0.580162492432517
	z = 0.4114188
	m = bubbleimg.imgmeasure.iso.isoMeasurer(dir_obj = 'verification_data_tabtools/SDSSJ0920+0034/', survey='hsc', z=z, center_mode='n/2-1')
	return m


def test_tabtools_has_row():

	fn_in = dir_test+'msr_iso.csv'
	assert tabtools.fn_has_row(fn_in, condi={'imgtag': 'OIII5008_I'})
	assert tabtools.fn_has_row(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})
	assert tabtools.fn_has_row(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'})


def test_tabtools_delete_row():
	fn_in = dir_test+'msr_iso.csv'
	assert tabtools.fn_has_row(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})

	tabtools.fn_delete_row(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})

	assert (not tabtools.fn_has_row(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'}))

	tab = at.Table.read(fn_in)
	assert len(tab) == 1


def test_tabtools_write_row_empty():

	fn_toadd = dir_verif+'msr_iso_toadd.csv'
	fn_in = dir_verif+'msr_iso.csv'
	fn_test = dir_test+'msr_iso_test.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	# write to an empty file
	assert not os.path.isfile(fn_test)
	tabtools.write_row(fn_test, row=tab_toadd, condi=condi, overwrite=False)
	assert os.path.isfile(fn_test)
	tab_test = at.Table.read(fn_test)
	assert len(tab_test) == 1 


def test_tabtools_write_row_nooverwrite():

	fn_toadd = dir_test+'msr_iso_toadd.csv'
	fn_in = dir_test+'msr_iso.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	assert len(at.Table.read(fn_in)) == 2
	# not overwrite an exisitng row
	tabtools.write_row(fn_in, row=tab_toadd, condi=condi, overwrite=False)
	assert os.path.isfile(fn_in)
	tab = at.Table.read(fn_in)
	assert tab['area_kpc'][1] > 0
	assert len(tab) == 2


def test_tabtools_write_row_overwrite():

	fn_toadd = dir_test+'msr_iso_toadd.csv'
	fn_in = dir_test+'msr_iso.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	assert len(at.Table.read(fn_in)) == 2
	# overwrite an exisitng row
	tabtools.write_row(fn_in, row=tab_toadd, condi=condi, overwrite=True)
	assert os.path.isfile(fn_in)
	tab = at.Table.read(fn_in)
	assert tab['area_kpc'][1] == 0
	assert len(tab) == 2


def test_tabtools_tab_extract_row():

	fn = dir_test+'msr_iso.csv'

	tab = at.Table.read(fn)

	assert len(tab) == 2

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}
	tab_ext = tabtools.tab_extract_row(tab, condi=condi)

	assert len(tab_ext) == 1



def test_isomeasurer_summarize(measurer1):
	m = measurer1
	fn = m.dir_obj+'msr_iso.csv'
	fn_sum = m.dir_obj+'msr_iso_smr.csv'

	imgtag = 'OIII5008_I'
	minarea = 5
	isocut1 = 1.e-15*u.Unit('erg / (arcsec2 cm2 s)')
	isocut2 = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')

	for isocut in (isocut1, isocut2):
		status = m.make_measurements(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=True, centerradius=5.*u.arcsec, overwrite=True, savecontours=False, plotmsr=False)

		assert status

	tab = at.Table.read(fn)
	assert len(tab) == 2

	condi = {'imgtag': imgtag}
	columns = ['area_ars']
	m.summarize(columns=columns, condi=condi, msrsuffix='')

	assert os.path.isfile(fn_sum)
	tab_sum = at.Table.read(fn_sum)
	assert len(tab_sum) == 1

	assert tab_sum['area_ars_mean'] == np.mean(tab['area_ars'])
	assert tab_sum['area_ars_std'] == np.std(tab['area_ars'])
	assert tab_sum['area_ars_median'] == np.median(tab['area_ars'])
	assert tab_sum['area_ars_p16'] == np.percentile(tab['area_ars'], 16)
	assert tab_sum['area_ars_p84'] == np.percentile(tab['area_ars'], 84)


def test_extract_line_from_file():

	fn_out = dir_test+'msr_iso_compile.csv'
	fn = dir_test +'hsc_xid.csv'
	fn_empty = dir_test+'hsc_xid_empty.csv'

	header = tabtools.extract_line_from_file(fn, iline=0)	
	l = tabtools.extract_line_from_file(fn, iline=1, comment='#', fill_trailing_empty=True)

	l_empty = tabtools.extract_line_from_file(fn_empty, iline=1, comment='#', fill_trailing_empty=True)

	tab_data = ascii.read([header, l, l_empty])
	tab_data.write(fn_out, format='ascii.csv')
	
	assert tab_data[1]['object_id'].mask == True
