# test_hscbatch_build.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at
import copy
import numpy as np

from ..hscbatch import hscBatch
from .... import imgdownload
from .... import obsobj



dir_parent = 'testing/'
dir_batch = 'testing/batch_ri/'
dir_batch_bad = 'testing/batch_ri_bad/'
dir_batch_exc = 'testing/batch_ri_exc/'
dir_batch_confus = 'testing/batch_ri_confus/'
name = 'batch_ri'
fn_cat = 'test_verification_data/example_catalog.fits'
fn_cat_bad = 'test_verification_data/bad_catalog.fits' # the last row of bad cat has bad ra, dec
fn_cat_exc = 'test_verification_data/only_exception_catalog.fits'
fn_cat_confus = 'test_verification_data/sandy_hsc_ra_confusion.fits'

fn_batch_photo_veri = 'test_verification_data/batch_hscphoto/' # the last row of bad cat has bad ra, dec
fn_batch_photo_test = 'testing/batch_hscphoto/' # the last row of bad cat has bad ra, dec


survey ='hsc'

fn_tab = 'sdss_photoobj.csv'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	shutil.copytree(fn_batch_photo_veri, fn_batch_photo_test)

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


@pytest.fixture
def batch1():
	return hscBatch(dir_batch=dir_batch, fn_cat=fn_cat)


@pytest.fixture
def batch_bad():
	return hscBatch(dir_batch=dir_batch_bad, fn_cat=fn_cat_bad)


@pytest.fixture
def batch_exc():
	return hscBatch(dir_batch=dir_batch_exc, fn_cat=fn_cat_exc)

@pytest.fixture
def batch_ra_confusion():
	return hscBatch(dir_batch=dir_batch_confus, fn_cat=fn_cat_confus)

@pytest.fixture
def batch_hscphotoobj_incomplete():
	return hscBatch(dir_batch=fn_batch_photo_test, fn_cat=fn_batch_photo_test+'list.csv', survey=survey)



def test_batch_compile_table(batch1):
	""" check compiled file exists and have correct content """
	b = batch1

	# build
	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status

	# compile table
	b.compile_table(fn_tab)

	fp = b.dir_batch+fn_tab
	assert os.path.isfile(fp)

	tab = at.Table.read(fp)
	assert len(tab) == len(b.list)

	# correct content
	for i in range(len(tab)):
		fp = b.dir_batch+'good/'+tab['obj_name'][i]+'/'+fn_tab
		objID = at.Table.read(fp)['objID'][0]

		assert tab['objID'][i] == objID


def test_batch_compile_table_w_except(batch_bad):
	b = batch_bad

	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status

	# compile table
	b.compile_table(fn_tab)

	# file exists
	fp = b.dir_batch+fn_tab
	assert os.path.isfile(fp)

	# file length correct
	tab = at.Table.read(fp)
	assert len(tab) == len(b.list)


	# correct content
	for i in range(len(tab)):
		if tab['obj_name'][i] in b.list_good['obj_name']: # if its a good object

			fp = b.dir_batch+'good/'+tab['obj_name'][i]+'/'+fn_tab
			objID = at.Table.read(fp)['objID'][0]

			assert tab['objID'][i] == objID
		else: 
			tab['obj_name'][i] == None


def test_batch_compile_table_w_only_except(batch_exc):
	b = batch_exc

	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status

	# compile table
	b.compile_table(fn_tab)

	fp = b.dir_batch+fn_tab
	assert not os.path.isfile(fp)



def test_batch_compile_table_ra_confusion(batch_ra_confusion):
	""" test whether the table compilation works when there are objects with very close RAs such taht the first part of the obj_names are identical """
	b = batch_ra_confusion

	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status

	# compile table
	b.compile_table(fn_tab)

	# file exists
	fp = b.dir_batch+fn_tab
	assert os.path.isfile(fp)

	# file length correct
	tab = at.Table.read(fp)
	assert len(tab) == len(b.list)

	# correct content
	for i in range(len(tab)):
		if tab['obj_name'][i] in b.list_good['obj_name']: # if its a good object

			fp = b.dir_batch+'good/'+tab['obj_name'][i]+'/'+fn_tab
			objID = at.Table.read(fp)['objID'][0]

			assert tab['objID'][i] == objID
		else: 
			tab['obj_name'][i] == None


def test_batch_compile_incomplete_table(batch_hscphotoobj_incomplete):

	b = batch_hscphotoobj_incomplete
	fn = b.dir_batch+'hsc_photoobj.csv'
	fn1 = b.dir_batch+'good/'+b.list[0]['obj_name']+'/hsc_photoobj.csv'

	tab1 = at.Table.read(fn1)

	assert len(b.list) > 0
	assert len(b.list_good) > 0

	b.compile_table('hsc_photoobj.csv')

	assert os.path.isfile(fn)

	tab = at.Table.read(fn)

	assert len(tab) == 2

	assert len(tab.colnames) == len(tab1.colnames) + 3



def func_build(obj, overwrite=False, **kwargs):
	"""
	Params
	------
	obj
	overwrite=False

	**kwargs:
		environment='iaa'

	Return
	------
	status
	"""

	# setting
	environment = kwargs.pop('environment', 'iaa')
	humvi_bands = 'riz'

	# running
	L = imgdownload.hscimgLoader(obj=obj, environment=environment)

	statuss = [
				L.status, 
				L.add_obj_sdss(), 
				]

	return all(statuss)


