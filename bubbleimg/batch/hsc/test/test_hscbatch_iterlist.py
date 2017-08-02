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

from setpaths import *
from fixture_built_batch import *


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
	

def test_batch_build_wexcept(batch_wexcept):
	b = batch_wexcept

	list_good = at.Table.read(b.dir_batch+'good/list_good.csv')
	list_except = at.Table.read(b.dir_batch+'except/list_except.csv')

	# to be constructed
	assert len(list_good) == 3
	assert len(b.list_good) == 3
	assert len(list_except) == 1
	assert len(b.list_except) == 1


def func_iterlist(obj, overwrite=False, **kwargs):

	fn = obj.dir_parent+fn_testing
	print fn
	with open(fn, 'a') as f:
		f.write(obj.name+'\n')

	return True