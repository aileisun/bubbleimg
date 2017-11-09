# test_hscbatch.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at
import copy

from ..hscbatch import hscBatch
from .... import obsobj

from .setpaths import *

obj_naming_sys = 'sdss'
catalog = at.Table.read(fn_cat, format='fits')


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
def batch_good():
	return hscBatch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)


def test_batch_init():
	"""
	init hscBatch with dir_batch (or dir_parent + name), catlogue, survey 
	"""

	b = hscBatch(dir_batch=dir_batch, catalog=catalog, survey=survey)
	assert b.name == name

	b = hscBatch(dir_parent=dir_parent, name=name, catalog=catalog, survey=survey)
	assert b.dir_batch == dir_batch

	assert b.survey == survey
	assert b.obj_naming_sys == obj_naming_sys
	assert len(b.catalog) == len(catalog)


def test_batch_init_error_noname():
	with pytest.raises(Exception):
		b = hscBatch(dir_parent=dir_parent, catalog=catalog, survey=survey)


def test_batch_init_error_badname():
	with pytest.raises(Exception):
		b = hscBatch(dir_batch=dir_batch, name='bad_name', catalog=catalog, survey=survey)


def test_batch_init_error_badcatalog():
	bad_catalog = copy.copy(catalog)
	bad_catalog.rename_column('RA', 'RA_new')	
	with pytest.raises(Exception):
		b = hscBatch(dir_batch=dir_batch, catalog=bad_catalog, survey=survey)


def test_batch_init_with_fn_cat():
	b = hscBatch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)
	assert len(b.catalog) == len(catalog)


def test_batch_write_list(batch_good):
	b = batch_good

	b._write_a_list(listname='')

	lst = at.Table.read(b.dir_batch+'list.csv')
	assert len(lst) > 0


def test_batch_args_to_list():

	args_to_list = ['z']

	b = hscBatch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey, args_to_list=args_to_list)

	b._write_a_list(listname='')
	lst = at.Table.read(b.dir_batch+'list.csv')

	for arg in args_to_list:
		assert arg in lst.colnames

def test_batch_override_short_obj_name_in_catalog():
	""" have obj_naming_sys overriding 'obj_name' in the input catalog """

	fn_list_w_short_obj_name = dir_verifc+'list_w_short_obj_name.csv'
	catalog = at.Table.read(fn_list_w_short_obj_name, format='ascii.csv')

	dir_batch = dir_parent + 'batch_short_obj_name/'
	survey = 'hsc'

	b = hscBatch(dir_batch=dir_batch, catalog=catalog, survey=survey)
	print((b.list[0]['obj_name']))
	assert len(b.list[0]['obj_name']) == 14

	b = hscBatch(dir_batch=dir_batch, catalog=catalog, survey=survey, obj_naming_sys='sdss_precise')

	print((b.list[0]['obj_name']))
	assert len(b.list[0]['obj_name']) == 18
