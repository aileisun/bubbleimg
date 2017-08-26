# test_batch.py
# ALS 2017/05/29


import astropy.table as at
import pytest
import os
import shutil
import copy
from astropy.io import ascii

from .. import batch
from ..batch import Batch
from ... import imgdownload
from ... import obsobj


dir_parent = 'testing/'
dir_batch = 'testing/batch_ri/'
dir_veri = 'test_verification_data/' 
name = 'batch_ri'
fn_cat = dir_veri+'example_catalog.fits'
fn_cat_bad = dir_veri+'bad_catalog.fits' # the last row of bad cat has bad ra, dec
catalog = at.Table.read(fn_cat, format='fits')
survey ='hsc'
obj_naming_sys = 'sdss'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def batch1():
	return Batch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)
	

@pytest.fixture
def batch_bad():
	return Batch(dir_batch=dir_batch, fn_cat_bad=fn_cat_bad, survey=survey)


def test_batch_init():
	"""
	init batch with dir_batch (or dir_parent + name), catlogue, survey 
	"""

	b = Batch(dir_batch=dir_batch, catalog=catalog, survey=survey)
	assert b.name == name

	b = Batch(dir_parent=dir_parent, name=name, catalog=catalog, survey=survey)
	assert b.dir_batch == dir_batch

	assert b.survey == survey
	assert b.obj_naming_sys == obj_naming_sys
	assert len(b.catalog) == len(catalog)



def test_batch_init_error_noname():
	with pytest.raises(Exception):
		b = Batch(dir_parent=dir_parent, catalog=catalog, survey=survey)


def test_batch_init_error_badname():
	with pytest.raises(Exception):
		b = Batch(dir_batch=dir_batch, name='bad_name', catalog=catalog, survey=survey)


def test_batch_init_error_badcatalog():
	bad_catalog = copy.copy(catalog)
	bad_catalog.rename_column('RA', 'RA_new')	
	with pytest.raises(Exception):
		b = Batch(dir_batch=dir_batch, catalog=bad_catalog, survey=survey)


def test_batch_init_with_fn_cat():
	b = Batch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)
	assert len(b.catalog) == len(catalog)


def test_batch_write_list(batch1):
	b = batch1

	b._write_a_list()

	lst = at.Table.read(b.dir_batch+'list.csv')
	assert len(lst) > 0


def test_extract_line_from_file():
	if not os.path.isdir(dir_parent):
		os.mkdir(dir_parent)

	fn_out = dir_parent+'hsc_xid_compiled.csv'
	fn = dir_veri +'hsc_xid.csv'
	header = batch._extract_line_from_file(fn, iline=0)	
	l = batch._extract_line_from_file(fn, iline=1, comment='#', fill_trailing_empty=True)

	fn_empty = dir_veri+'hsc_xid_empty.csv'
	l_empty = batch._extract_line_from_file(fn_empty, iline=1, comment='#', fill_trailing_empty=True)

	tab_data = ascii.read([header, l, l_empty])
	tab_data.write(fn_out, format='ascii.csv')
	
	assert tab_data[1]['object_id'].mask == True
