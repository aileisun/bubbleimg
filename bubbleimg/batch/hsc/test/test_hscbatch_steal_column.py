# test_hscbatch_steal_columns.py
# ALS 2017/05/29

import pytest
import os
import shutil
import astropy.table as at

from setpaths import *
from fixture_built_batch import *

fn_tab = 'sdss_xid.csv'


@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)
	os.mkdir(dir_parent)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


def test_batch_steal_columns_good(batch_good):
	""" check compiled file exists and have correct content """

	b = batch_good
	col='patch_id'

	tab_cat = at.Table.read(fn_cat)
	tab_cat.rename_column('RA', 'ra')
	tab_cat.rename_column('DEC', 'dec')
	tab_cat = tab_cat[::-1]
	assert col in tab_cat.colnames

	b.steal_columns(tab=tab_cat, colnames=[col])

	for lst in [b.list, b.list_good, b.list_except, ]:
		assert col in lst.colnames

		if len(lst)>0:
			print lst
			print tab_cat
			tab_join = at.join(lst, tab_cat, keys=['ra', 'dec'], join_type='left')
			assert all(tab_join[col+'_1'] == tab_join[col+'_2'])


	for fn in [b.fp_list, b.fp_list_good, b.fp_list_except]:
		tab = at.Table.read(fn)
		assert col in tab.colnames


