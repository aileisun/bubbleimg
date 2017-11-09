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

from .setpaths import *
from .fixture_built_batch import *


fn_testing = 'testing_list.txt'

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


def test_batch_iterlist(batch_good):
	b = batch_good

	fn = dir_batch+'good/'+fn_testing
	if os.path.isfile(fn):
		os.remove(fn)

	statuss = b.iterlist(func_iterlist, listname='good')

	assert all(statuss)

	assert os.path.isfile(fn)

	a = np.genfromtxt(fn)

	assert a.sort() == b.list_good['obj_name'].sort()
	

def test_batch_iterlist_imgdownload_squential(batch_good):
	b = batch_good

	statuss = b.iterlist(func_iterlist_imgdownload, listname='good', overwrite=True, processes=-1)

	assert all(statuss)

	obj = b.get_ith_obj_from_list(iobj=0, listname='good')

	assert os.path.isfile(obj.dir_obj + 'stamp-i.fits')



def test_batch_iterlist_imgdownload_3processes(batch_good):
	b = batch_good

	statuss = b.iterlist(func_iterlist_imgdownload, listname='good', overwrite=True, processes=3)

	assert all(statuss)

	obj = b.get_ith_obj_from_list(iobj=0, listname='good')

	assert os.path.isfile(obj.dir_obj + 'stamp-i.fits')


def test_batch_iterlist_imgdownload_default_processes(batch_good):
	b = batch_good

	statuss = b.iterlist(func_iterlist_imgdownload, listname='good', overwrite=True, processes=None)

	assert all(statuss)

	obj = b.get_ith_obj_from_list(iobj=0, listname='good')

	assert os.path.isfile(obj.dir_obj + 'stamp-i.fits')


def func_iterlist(obj, overwrite=False, **kwargs):

	fn = obj.dir_parent+fn_testing
	print(fn)
	with open(fn, 'a') as f:
		f.write(obj.name+'\n')

	return True


def func_iterlist_imgdownload(obj, overwrite=False, **kwargs):

	obj.add_hsc(overwrite = overwrite)
	L = imgdownload.hscimgLoader(obj=obj)

	status = L.make_stamps(overwrite=overwrite)

	return status