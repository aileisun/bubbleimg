import pytest
import os
import astropy.table as at
import shutil

from ..sdssobj import sdssObj

ra = 150.0547735
dec = 12.7073027

dir_obj = './testing/SDSSJ1000+1242/'
dir_parent = './testing/'

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
	return sdssObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_SDSSObj_init_dir_obj(obj_dirobj):

	obj = obj_dirobj

	assert isinstance(obj, sdssObj)
	assert obj.dir_obj == dir_obj


def test_SDSSObj_init_dir_parent():

	obj = sdssObj(ra=ra, dec=dec, dir_parent=dir_parent)

	assert isinstance(obj, sdssObj)
	assert obj.dir_obj == dir_obj


def test_SDSSObj_xid(obj_dirobj):

	obj = obj_dirobj
	status = obj.load_xid(writefile=False)

	assert status
	assert round(obj.xid['ra'], 4) == round(ra, 4)


def test_SDSSObj_xid_writefile(obj_dirobj):

	obj = obj_dirobj
	fn = obj.dir_obj+'sdss_xid.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_xid(writefile=True)

	assert status
	assert round(obj.xid['ra'], 4) == round(ra, 4)

	assert os.path.isfile(fn)

	fxid = at.Table.read(fn, format='csv')
	assert round(fxid['ra'][0], 4) == round(ra, 4)


def test_SDSSObj_xid_fails():
	ra = 0.
	dec = -89.
	obj = sdssObj(ra=ra, dec=dec, dir_obj = './badobject/')

	status = obj.load_xid(writefile=True)

	assert status == False

	fn = obj.dir_obj+'sdss_xid.csv'
	assert not os.path.isfile(fn)


def test_SDSSObj_xid_conflicting_dir_obj():
	# use a dir_obj that is occupied by another different obj
	ra = 150.0547735
	dec = 12.7073027

	dir_obj = './testing/SDSSJ9999+9999/'

	with pytest.raises(Exception):
		obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj)


def test_load_photoboj_writefile(obj_dirobj):
	obj = obj_dirobj
	fn = obj.dir_obj+'sdss_photoobj.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_photoobj(writefile=True)
	assert os.path.isfile(fn)
	assert status


def test_SDSSObj_photoobj_fails():
	ra = 0.
	dec = -89.
	obj = sdssObj(ra=ra, dec=dec, dir_obj = './badobject/')

	status = obj.load_photoobj(writefile=True)

	assert status == False

	fn = obj.dir_obj+'sdss_photoobj.csv'
	assert not os.path.isfile(fn)




def test_SDSSObj_make_spec_writefile(obj_dirobj):
	obj = obj_dirobj
	fn = obj.dir_obj+'spec.fits'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	spec = obj.make_spec(overwrite=True)
	assert os.path.isfile(fn)


def test_SDSSObj_get_speclcoord(obj_dirobj):
	obj = obj_dirobj
	spec, lcoord = obj.get_speclcoord(wunit=False)

	assert len(spec) > 0
	assert len(spec) == len(lcoord)
	assert sum(spec) > 0


def test_SDSSObj_get_spec(obj_dirobj):
	obj = obj_dirobj
	spechdu = obj.get_spec()

	assert len(spechdu[1].data['flux']) > 1