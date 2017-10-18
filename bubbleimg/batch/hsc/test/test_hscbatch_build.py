# test_hscbatch_build.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at
import numpy as np
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

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def batch_good():
	return hscBatch(dir_batch=dir_batch, fn_cat=fn_cat)


@pytest.fixture
def batch_wexcept():
	return hscBatch(dir_batch=dir_batch_wexcept, fn_cat=fn_cat_wexcept)


def test_batch_build(batch_good):
	b = batch_good

	kwargs = {'environment': 'online'}
	status = b.build(func_build, processes=-1, **kwargs)

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


def test_batch_build_other_batches_sequential():

	b1 = hscBatch(dir_batch=dir_batch_onlyexcept, fn_cat=fn_cat_onlyexcept)
	b2 = hscBatch(dir_batch=dir_batch_confus, fn_cat=fn_cat_confus)

	for b in [b1, b2]:
		status = b.build(func_build_hsc_sdss, processes=-1)
		assert status

		check_lst_in_order(b.list)
		check_lst_in_order(b.list_good)
		check_lst_in_order(b.list_except)


def test_batch_build_other_batches_default_processes():

	b1 = hscBatch(dir_batch=dir_batch_onlyexcept, fn_cat=fn_cat_onlyexcept)
	b2 = hscBatch(dir_batch=dir_batch_confus, fn_cat=fn_cat_confus)

	for b in [b1, b2]:
		status = b.build(func_build_hsc_sdss, processes=None)
		assert status

		check_lst_in_order(b.list)
		check_lst_in_order(b.list_good)
		check_lst_in_order(b.list_except)


def test_batch_build_except_twice():
	b = hscBatch(dir_batch=dir_batch_wexcept, fn_cat=fn_cat_wexcept)

	status = b.build(func_build_hsc_sdss, overwrite=True)
	assert status
	assert os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/')
	assert not os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/SDSSJ0000-8900/')

	status = b.build(func_build_hsc_sdss, overwrite=False)
	assert status
	assert os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/')
	assert not os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/SDSSJ0000-8900/')

	status = b.build(func_build_hsc_sdss, overwrite=True)
	assert status
	assert os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/')
	assert not os.path.isdir(b.dir_batch+'except/SDSSJ0000-8900/SDSSJ0000-8900/')


def check_lst_in_order(lst):

	lst_sort = copy.deepcopy(lst)
	lst_sort.sort('ra')

	assert np.all(lst == lst_sort)


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
	statuss = [obj.add_hsc(), 
				obj.add_sdss()
				]

	return all(statuss)


