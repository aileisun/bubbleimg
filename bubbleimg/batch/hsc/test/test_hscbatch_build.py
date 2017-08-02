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

from setpaths import *


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
def batch_good():
	return hscBatch(dir_batch=dir_batch, fn_cat=fn_cat)


@pytest.fixture
def batch_wexcept():
	return hscBatch(dir_batch=dir_batch_wexcept, fn_cat=fn_cat_wexcept)


def test_batch_build(batch_good):
	b = batch_good

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


def test_batch_build_wexcept(batch_wexcept):
	b = batch_wexcept

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


def test_batch_build_check_folders_consistent_w_list(batch_good):
	b = batch_good

	status = b.build(func_build_hscobj)
	b._check_folders_consistent_w_list()

	dir_obj = b.dir_good+b.list_good['obj_name'][-1]+'/'

	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)

	with pytest.raises(Exception):
		b._check_folders_consistent_w_list()	

def test_batch_build_other_batches():

	b1 = hscBatch(dir_batch=dir_batch_onlyexcept, fn_cat=fn_cat_onlyexcept)
	b2 = hscBatch(dir_batch=dir_batch_confus, fn_cat=fn_cat_confus)

	for b in [b1, b2]:
		status = b.build(func_build_hsc_sdss)
		assert status


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


def func_build_hsc_sdss(obj, overwrite=False, **kwargs):
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


