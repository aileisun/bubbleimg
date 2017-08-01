# test_hscbatch_build.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at
import copy

from ..hscbatch import hscBatch
from .... import imgdownload
from .... import obsobj


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





def test_batch_build_check_folders_consistent_w_list(batch1):
	b = batch1

	status = b.build(func_build_hscobj)
	b._check_folders_consistent_w_list()

	dir_obj = b.dir_good+b.list_good['obj_name'][-1]+'/'

	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)

	with pytest.raises(Exception):
		b._check_folders_consistent_w_list()

	


def func_build_hscobj(obj, overwrite=False):
	"""
	Params
	------
	obj	
	overwrite=False

	Return
	------
	status
	"""

	return obj.add_hsc()


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
	humvi_bands = 'riz'

	# running
	L = imgdownload.hscimgLoader(obj=obj, **kwargs)

	if L.status:
		statuss = 	[ 
					L.make_stamps(overwrite=overwrite), 
					L.make_psfs(overwrite=overwrite), 
					L.plot_colorimg(bands=humvi_bands, img_type='stamp', overwrite=overwrite),
					L.add_obj_sdss(), 
					L.obj.sdss.make_spec(overwrite=overwrite),
					]

		return all(statuss)
	else:
		return False
