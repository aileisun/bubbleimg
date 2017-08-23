import pytest
import os
import astropy.table as at
import numpy as np
import shutil
import filecmp
from astropy.io import ascii

from ..hscobj import hscObj

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

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


@pytest.fixture
def obj_dirobj():
	if os.path.isdir(dir_obj):
		shutil.rmtree(dir_obj)
	return hscObj(ra=ra, dec=dec, dir_obj=dir_obj, rerun='s16a_wide2')


@pytest.fixture
def obj_dirobj_twohscsource():
	ra = 18.3349643407323022
	dec = 2.85719938197981449
	return hscObj(ra=ra, dec=dec, dir_parent=dir_parent, rerun='s16a_wide2')


def test_hscObj_init_dir_obj(obj_dirobj):

	obj = obj_dirobj

	assert isinstance(obj, hscObj)
	assert obj.dir_obj == dir_obj


def test_hscObj_init_dir_parent():

	obj = hscObj(ra=ra, dec=dec, dir_parent=dir_parent)

	assert isinstance(obj, hscObj)
	assert obj.dir_obj == dir_obj


def test_hscObj_xid(obj_dirobj):

	obj = obj_dirobj
	status = obj.load_xid()

	assert status
	assert np.absolute(obj.xid['ra'][0] - ra) < 0.0003
	assert 'tract' in obj.xid.colnames
	assert 'patch' in obj.xid.colnames
	assert hasattr(obj, 'tract')
	assert hasattr(obj, 'patch')

	assert obj.tract == 9564
	assert obj.patch == 703

	fn = obj.dir_obj+'hsc_xid.csv'
	fxid = at.Table.read(fn, format='csv')
	assert np.absolute(fxid['ra'][0] - ra) < 0.0003
	assert fxid['tract'][0] == 9564


def test_hscObj_xid_overwrite(obj_dirobj):

	obj = obj_dirobj
	fn = obj.dir_obj+'hsc_xid.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_xid(overwrite=True)

	assert os.path.isfile(fn)
	
	# when file does exist it overwrites stamp
	if os.path.isfile(fn):
		os.remove(fn)
	open(fn, 'w').close()
	assert os.stat(fn).st_size == 0

	with pytest.raises(Exception):
		status = obj.load_xid(overwrite=False)
	assert os.stat(fn).st_size == 0

	status = obj.load_xid(overwrite=True)
	assert os.path.isfile(fn)
	assert os.stat(fn).st_size > 0


def test_hscObj_xid_fails():
	ra = 0.
	dec = -89.
	obj = hscObj(ra=ra, dec=dec, dir_obj = './testing/badobject/')

	status = obj.load_xid()

	assert status == False

	fn = obj.dir_obj+'sdss_xid.csv'
	assert not os.path.isfile(fn)


def test_hscObj_xid_conflicting_dir_obj():
	# use a dir_obj that is occupied by another different obj
	ra = 140.099341430207
	dec = 0.580162492432517

	dir_obj = './testing/SDSSJ9999+9999/'

	with pytest.raises(Exception):
		obj = hscObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_only_one_row(obj_dirobj):
	obj = obj_dirobj

	assert len(obj.xid) == 1

	tab = at.Table.read(obj.dir_obj+'hsc_xid.csv', format='ascii.csv')
	assert len(tab) == 1


def test_only_one_row_twohscsource(obj_dirobj_twohscsource):
	obj = obj_dirobj_twohscsource

	assert len(obj.xid) == 1

	tab = at.Table.read(obj.dir_obj+'hsc_xid.csv', format='ascii.csv')
	assert len(tab) == 1


def test_hscObj_identical_w_verification(obj_dirobj):

	f = 'hsc_xid.csv'
	file_totest = './testing/SDSSJ0920+0034/'+f
	file_verification = './test_verification_data/SDSSJ0920+0034/'+f

	assert filecmp.cmp(file_totest, file_verification)


def test_hscObj_download_photoobj(obj_dirobj):
	obj = obj_dirobj
	assert obj.status

	band_columns = ['mag_kron', 'mag_kron_err', 'flux_kron_flags', 'flux_kron_radius', 'mag_aperture10', 'mag_aperture15']
	bands = ['g', 'r', 'i', 'z', 'y'] 

	status = obj._download_photoobj(band_columns=band_columns, bands=bands, all_columns=False, tab_name='forced', overwrite=True)

	assert status
	photoobj = at.Table.read(obj.dir_obj+'hsc_photoobj.csv')

	assert len(photoobj) == 1
	assert len(photoobj.colnames) > 1

	for x in bands:
		for y in band_columns:
			assert x+y in photoobj.colnames
	

def test_hscObj_loadphotoobj(obj_dirobj):

	band_columns = ['mag_kron', 'mag_kron_err', 'flux_kron_flags', 'flux_kron_radius', 'mag_aperture10', 'mag_aperture15']
	bands = ['g', 'r', 'i', 'z', 'y'] 


	obj = obj_dirobj
	assert obj.status

	status = obj.load_photoobj(band_columns=[], bands=[], tab_name='forced', all_columns=False, overwrite=True)

	assert status
	assert os.path.isfile(obj.fp_photoobj)
	assert len(obj.photoobj) == 1
	assert len(obj.photoobj.colnames) > 1

	for col in band_columns:
		for b in bands:
			assert b+col in obj.photoobj.colnames

	band_columns = ['mag_kron', 'mag_kron_err']
	status = obj.load_photoobj(band_columns=band_columns, bands=[], tab_name='forced', all_columns=False, overwrite=True)

	assert status
	assert os.path.isfile(obj.fp_photoobj)
	assert len(obj.photoobj) == 1
	assert len(obj.photoobj.colnames) > 1

	for col in band_columns:
		for b in bands:
			assert b+col in obj.photoobj.colnames


def test_hscObj_loadphotoobj_catalog(obj_dirobj):
	obj = obj_dirobj
	assert obj.status
	status = obj.load_photoobj(band_columns=[], bands=[], tab_name='forced', all_columns=False, overwrite=True)
	assert status
	assert obj.photoobj['tab_name'] == 'forced'


	status = obj.load_photoobj(band_columns=[], bands=[], tab_name='testing', all_columns=False, overwrite=True)
	assert status == False


def test_hscObj_xid_sanity_check(obj_dirobj_twohscsource):
	obj = obj_dirobj_twohscsource

	tabstr = ['object_id,ra,dec,patch_id,tract,patch,patch_s,parent_id,deblend_nchild,detect_is_patch_inner,detect_is_tract_inner,detect_is_primary', '42766771777714198,18.3349643407323022,2.85719938197981449,97240107,9724,107,"1,7",42766771777704241,0,t,t,t', '42766771777714199,18.3353134708619052,2.85678842440856107,97240107,9724,107,"1,7",42766771777704241,0,t,t,t',]

	xid_tworows = ascii.read(tabstr)
	xid = at.Table(xid_tworows[0])

	# more than one instance
	with pytest.raises(Exception):
		obj._xid_sanity_check(xid_tworows)

	# zero instance
	with pytest.raises(Exception):
		obj._xid_sanity_check(at.Table())

	# ok xid
	obj._xid_sanity_check(xid)

	# is not primary
	xid['detect_is_primary'][0] = 'f'
	with pytest.raises(Exception):
		obj._xid_sanity_check(xid)

	# ra dec wrong
	xid['detect_is_primary'][0] = 't'
	xid['ra'][0] = xid['ra'][0] + 0.05
	with pytest.raises(Exception):
		obj._xid_sanity_check(xid)


def test_hscObj_photoobj_overwrite(obj_dirobj):

	obj = obj_dirobj
	fn = obj.dir_obj+'hsc_photoobj.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_photoobj(overwrite=True)

	assert os.path.isfile(fn)
	assert obj.photoobj['object_id'][0] == 42063891789801134

	tab = at.Table.read(fn)
	tab['object_id'][0] = 0
	tab.write(fn, format='ascii.csv', overwrite=True)
	
	status = obj.load_photoobj(overwrite=False)
	assert obj.photoobj['object_id'][0] == 0
	tab = at.Table.read(fn)
	assert tab['object_id'][0] == 0

	status = obj.load_photoobj(overwrite=True)
	assert obj.photoobj['object_id'][0] == 42063891789801134
	tab = at.Table.read(fn)
	assert tab['object_id'][0] == 42063891789801134


def test_hscObj_download_table(obj_dirobj):
	obj = obj_dirobj
	fn = obj.dir_obj+'hsc_photoz_demp.csv'

	status = obj.download_table(tab_name='photoz_demp', columns=[], overwrite=True)
	assert status

	assert os.path.isfile(fn)

	tab = at.Table.read(fn, format='ascii.csv')

	assert 'photoz_best' in tab.colnames
	assert len(tab) == 1

