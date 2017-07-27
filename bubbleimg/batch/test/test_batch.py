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

fn_batch_photo_veri = 'test_verification_data/batch_hscphoto/' # the last row of bad cat has bad ra, dec
fn_batch_photo_test = 'testing/batch_hscphoto/' # the last row of bad cat has bad ra, dec


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	shutil.copytree(fn_batch_photo_veri, fn_batch_photo_test)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def batch1():
	return Batch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)


@pytest.fixture
def batch_hscphotoobj():
	return Batch(dir_batch=fn_batch_photo_test, fn_cat=fn_batch_photo_test+'list.csv', survey=survey)


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


def test_batch_compile_table(batch_hscphotoobj):

	b = batch_hscphotoobj
	fn = b.dir_batch+'hsc_photoobj.csv'

	assert len(b.list) > 0
	assert len(b.list_good) > 0

	b.compile_table('hsc_photoobj.csv')

	assert os.path.isfile(fn)

	tab = at.Table.read(fn)

	assert len(tab) == 2

