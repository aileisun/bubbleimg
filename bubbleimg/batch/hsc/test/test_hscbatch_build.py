# test_hscbatch_build.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at
import copy

from ..hscbatch import hscBatch
from .... import downloadimg
from .... import obsobj
from .... import blobmaps


dir_parent = 'testing/'
dir_batch = 'testing/batch_ri/'
dir_batch_bad = 'testing/batch_ri_bad/'
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
	return hscBatch(dir_batch=dir_batch, fn_cat=fn_cat)


@pytest.fixture
def batch_bad():
	return hscBatch(dir_batch=dir_batch_bad, fn_cat=fn_cat_bad)


def test_batch_build(batch1):
	b = batch1

	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status
	assert os.path.isdir(b.dir_batch+'good/')
	assert os.path.isdir(b.dir_batch+'except/')

	assert os.path.isfile(b.dir_batch+'good/list_good.csv')
	assert os.path.isfile(b.dir_batch+'except/list_except.csv')

	list_good = at.Table.read(b.dir_batch+'good/list_good.csv')
	list_except = at.Table.read(b.dir_batch+'except/list_except.csv')

	assert 'obj_name' in list_good.colnames
	assert len(list_good) == 3
	assert len(list_except) == 0

	for obj_name in list_good['obj_name']:
		assert os.path.isdir(b.dir_batch+'good/'+obj_name+'/')

	for obj_name in list_except['obj_name']:
		assert os.path.isdir(b.dir_batch+'except/'+obj_name+'/')


def test_batch_build_bad(batch_bad):
	b = batch_bad

	kwargs = {'environment': 'online'}
	status = b.build(func_build, **kwargs)

	assert status
	assert os.path.isdir(b.dir_batch+'good/')
	assert os.path.isdir(b.dir_batch+'except/')

	assert os.path.isfile(b.dir_batch+'good/list_good.csv')
	assert os.path.isfile(b.dir_batch+'except/list_except.csv')

	list_good = at.Table.read(b.dir_batch+'good/list_good.csv')
	list_except = at.Table.read(b.dir_batch+'except/list_except.csv')

	assert 'obj_name' in list_good.colnames
	assert len(list_good) == 3
	assert len(list_except) == 1

	for obj_name in list_good['obj_name']:
		assert os.path.isdir(b.dir_batch+'good/'+obj_name+'/')

	for obj_name in list_except['obj_name']:
		assert os.path.isdir(b.dir_batch+'except/'+obj_name+'/')



def func_build(obj, overwrite=False, **kwargs):
	"""
	Params
	------
	ra
	dec
	obj_name
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
	L = downloadimg.hscimgLoader(obj=obj, environment=environment)

	statuss = [
				L.hsc_status, 
				L.add_obj_sdss(), 
				L.make_stamps(overwrite=overwrite), 
				L.plot_colorimg(bands=humvi_bands, img_type='stamp', overwrite=overwrite), 
				]

	return all(statuss)