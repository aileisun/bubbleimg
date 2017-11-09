# test_hscbatch_compile_table.py
# ALS 2017/05/29


import pytest
import os
import shutil
import astropy.table as at

from .setpaths import *
from .fixture_built_batch import *

fn_tab = 'sdss_xid.csv'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)
	os.mkdir(dir_parent)

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


def test_batch_compile_table(batch_good):
	""" check compiled file exists and have correct content """
	b = batch_good

	# compile table
	b.compile_table(fn_tab)

	fp = b.dir_batch+fn_tab
	assert os.path.isfile(fp)

	tab = at.Table.read(fp)
	assert len(tab) == len(b.list)

	# correct content
	for i in range(len(tab)):
		fp = b.dir_batch+'good/'+tab['obj_name'][i]+'/'+fn_tab
		objid = at.Table.read(fp)['objid'][0]

		assert tab['objid'][i] == objid


def test_batch_compile_table_w_except(batch_wexcept):
	b = batch_wexcept

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
			objid = at.Table.read(fp)['objid'][0]

			assert tab['objid'][i] == objid
		else: 
			tab['obj_name'][i] == None


def test_batch_compile_table_w_only_except(batch_onlyexcept):
	b = batch_onlyexcept

	# compile table
	b.compile_table(fn_tab)

	fp = b.dir_batch+fn_tab
	assert not os.path.isfile(fp)



def test_batch_compile_table_ra_confusion(batch_confus):
	""" test whether the table compilation works when there are objects with very close RAs such taht the first part of the obj_names are identical """
	b = batch_confus

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
			objid = at.Table.read(fp)['objid'][0]

			assert tab['objid'][i] == objid
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



def test_batch_compile_multiline(batch_multiline):

	b = batch_multiline
	fn = 'psf_smeared.csv'
	fp = b.dir_batch+fn

	b.compile_table(fn, alllines=True)

	assert os.path.isfile(fp)

	tab = at.Table.read(fp)

	assert len(tab) == 15

	assert 'psf_fwhm_arcs' in tab.colnames

