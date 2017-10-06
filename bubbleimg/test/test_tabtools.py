
import pytest
import os
import shutil
import astropy.table as at
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


def test_tabtools_has_line():

	fn_in = dir_test+'msr_iso.csv'
	assert tabtools.fn_has_line(fn_in, condi={'imgtag': 'OIII5008_I'})
	assert tabtools.fn_has_line(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})
	assert tabtools.fn_has_line(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'})


def test_tabtools_delete_line():
	fn_in = dir_test+'msr_iso.csv'
	assert tabtools.fn_has_line(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})

	tabtools.fn_delete_line(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'})

	assert (not tabtools.fn_has_line(fn_in, condi={'imgtag': 'OIII5008_I', 'isocut':'1e-15 erg / (arcsec2 cm2 s)'}))

	tab = at.Table.read(fn_in)
	assert len(tab) == 1


def test_tabtools_write_line_empty():

	fn_toadd = dir_test+'msr_iso_toadd.csv'
	fn_in = dir_test+'msr_iso.csv'
	fn_test = dir_test+'msr_iso_test.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	# write to an empty file
	assert not os.path.isfile(fn_test)
	tabtools.write_line(fn_test, line=tab_toadd, condi=condi, overwrite=False)
	assert os.path.isfile(fn_test)
	tab_test = at.Table.read(fn_test)
	assert len(tab_test) == 1 


def test_tabtools_write_line_nooverwrite():

	fn_toadd = dir_test+'msr_iso_toadd.csv'
	fn_in = dir_test+'msr_iso.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	assert len(at.Table.read(fn_in)) == 2
	# not overwrite an exisitng row
	tabtools.write_line(fn_in, line=tab_toadd, condi=condi, overwrite=False)
	assert os.path.isfile(fn_in)
	tab = at.Table.read(fn_in)
	assert tab['area_kpc'][1] > 0
	assert len(tab) == 2


def test_tabtools_write_line_overwrite(setUp_tearDown):

	fn_toadd = dir_test+'msr_iso_toadd.csv'
	fn_in = dir_test+'msr_iso.csv'

	tab_toadd = at.Table.read(fn_toadd)

	condi = {'imgtag': 'OIII5008_I', 'isocut':'3e-15 erg / (arcsec2 cm2 s)'}

	assert len(at.Table.read(fn_in)) == 2
	# overwrite an exisitng row
	tabtools.write_line(fn_in, line=tab_toadd, condi=condi, overwrite=True)
	assert os.path.isfile(fn_in)
	tab = at.Table.read(fn_in)
	assert tab['area_kpc'][1] == 0
	assert len(tab) == 2

