# test_batch.py
# ALS 2017/05/29


import astropy.table as at
import pytest
import os
import shutil
import copy

from ..batch import Batch
from ... import imgdownload
from ... import obsobj


dir_parent = 'testing/'
dir_batch = 'testing/batch_ri/'
name = 'batch_ri'
fn_cat = 'test_verification_data/example_catalog.fits'
fn_cat_bad = 'test_verification_data/bad_catalog.fits' # the last row of bad cat has bad ra, dec
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

	b._write_list()

	lst = at.Table.read(b.dir_batch+'list.csv')
	assert len(lst) > 0

