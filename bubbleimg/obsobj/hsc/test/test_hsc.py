import pytest
import os
import astropy.table as at
import numpy as np
import shutil
import filecmp

from ..hscobj import HSCObj

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'


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
def obj_dirobj():
	return HSCObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_HSCObj_init_dir_obj(obj_dirobj):

	obj = obj_dirobj

	assert isinstance(obj, HSCObj)
	assert obj.dir_obj == dir_obj


def test_HSCObj_init_dir_parent():

	obj = HSCObj(ra=ra, dec=dec, dir_parent=dir_parent)

	assert isinstance(obj, HSCObj)
	assert obj.dir_obj == dir_obj


def test_HSCObj_xid(obj_dirobj):

	obj = obj_dirobj
	status = obj.load_xid(writefile=False)

	assert status
	assert np.absolute(obj.xid['ra'][0] - ra) < 0.0003
	assert 'tract' in obj.xid.colnames
	assert 'patch' in obj.xid.colnames
	assert hasattr(obj, 'tract')
	assert hasattr(obj, 'patch')

	assert obj.tract == 9564
	assert obj.patch == 703


def test_HSCObj_xid_writefile(obj_dirobj):

	obj = obj_dirobj
	fn = obj.dir_obj+'hsc_xid.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_xid(writefile=True)

	assert status
	assert np.absolute(obj.xid['ra'][0] - ra) < 0.0003

	assert os.path.isfile(fn)

	fxid = at.Table.read(fn, format='csv')
	assert np.absolute(fxid['ra'][0] - ra) < 0.0003


def test_HSCObj_xid_fails():
	ra = 0.
	dec = -89.
	obj = HSCObj(ra=ra, dec=dec, dir_obj = './testing/badobject/')

	status = obj.load_xid(writefile=True)

	assert status == False

	fn = obj.dir_obj+'sdss_xid.csv'
	assert not os.path.isfile(fn)


def test_HSCObj_xid_conflicting_dir_obj():
	# use a dir_obj that is occupied by another different obj
	ra = 140.099341430207
	dec = 0.580162492432517

	dir_obj = './testing/SDSSJ9999+9999/'

	with pytest.raises(Exception):
		obj = HSCObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_only_one_row(obj_dirobj):
	obj = obj_dirobj

	assert len(obj.xid) == 1


def test_HSCObj_identical_w_verification(obj_dirobj):

	f = 'hsc_xid.csv'
	file_totest = './testing/SDSSJ0920+0034/'+f
	file_verification = './test_verification_data/SDSSJ0920+0034/'+f

	assert filecmp.cmp(file_totest, file_verification)
