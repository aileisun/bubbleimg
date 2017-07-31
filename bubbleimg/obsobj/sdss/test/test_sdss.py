import pytest
import os
import astropy.table as at
import shutil
import astropy.units as u

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

	# yield
	# # tear down
	# if os.path.isdir(dir_parent):
	# 	shutil.rmtree(dir_parent)


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


# def test_SDSSObj_xid(obj_dirobj):

# 	obj = obj_dirobj
# 	status = obj.load_xid(overwrite=True)

# 	assert status
# 	assert round(obj.xid['ra'], 4) == round(ra, 4)


def test_SDSSObj_xid(obj_dirobj):

	obj = obj_dirobj
	fn = obj.dir_obj+'sdss_xid.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_xid(overwrite=True)

	assert status
	assert round(obj.xid['ra'], 4) == round(ra, 4)

	assert os.path.isfile(fn)

	fxid = at.Table.read(fn, format='csv')
	assert round(fxid['ra'][0], 4) == round(ra, 4)


def test_SDSSObj_xid_overwrite():
	ra = 34.7489566190744
	dec = -4.04365821949229
	dir_obj = './testing/SDSSJ0218-0402/'

	run2d_dr12 = 'v5_7_0'
	run2d_dr13 = 'v5_9_0'

	if os.path.isfile(dir_obj):
		os.remove(dir_obj)
	assert not os.path.isfile(dir_obj)


	obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj, data_release=12, overwrite=True)

	fn = obj.dir_obj+'sdss_xid.csv'

	assert obj.status
	assert round(obj.xid['ra'], 4) == round(ra, 4)
	assert obj.run2d == run2d_dr12

	assert os.path.isfile(fn)

	fxid = at.Table.read(fn, format='csv')
	assert fxid['run2d'][0] == run2d_dr12


	obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj, data_release=13, overwrite=False)
	assert obj.run2d == run2d_dr12
	fxid = at.Table.read(fn, format='csv')
	assert fxid['run2d'][0] == run2d_dr12


	obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj, data_release=13, overwrite=True)

	assert obj.run2d == run2d_dr13
	fxid = at.Table.read(fn, format='csv')
	assert fxid['run2d'][0] == run2d_dr13


def test_SDSSObj_xid_fails():
	ra = 0.
	dec = -89.
	obj = sdssObj(ra=ra, dec=dec, dir_obj = './badobject/')

	status = obj.load_xid(overwrite=True)

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


def test_load_photoboj_overwrite(obj_dirobj):
	obj = obj_dirobj
	fn = obj.dir_obj+'sdss_photoobj.csv'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	status = obj.load_photoobj(overwrite=True)
	assert os.path.isfile(fn)
	assert status

	assert obj.photoobj['objID'][0] == 1237664106315579464

	tab = at.Table.read(fn)
	tab['objID'] = 0000000000
	tab.write(fn, overwrite=True)

	status = obj.load_photoobj(overwrite=False)
	assert obj.photoobj['objID'][0] == 0000000000

	status = obj.load_photoobj(overwrite=True)
	assert obj.photoobj['objID'][0] == 1237664106315579464


def test_SDSSObj_photoobj_fails():
	ra = 0.
	dec = -89.
	obj = sdssObj(ra=ra, dec=dec, dir_obj = './badobject/')

	status = obj.load_photoobj(overwrite=True)

	assert status == False

	fn = obj.dir_obj+'sdss_photoobj.csv'
	assert not os.path.isfile(fn)


def test_SDSSObj_make_spec_overwrite(obj_dirobj):
	obj = obj_dirobj
	fn = obj.dir_obj+'spec.fits'

	if os.path.isfile(fn):
		os.remove(fn)
	assert not os.path.isfile(fn)

	spec = obj.make_spec(overwrite=True)
	assert os.path.isfile(fn)

	# when file does exist it overwrites stamp
	if os.path.isfile(fn):
		os.remove(fn)
	open(fn, 'w').close()
	assert os.stat(fn).st_size == 0

	spec = obj.make_spec(overwrite=False)
	assert os.stat(fn).st_size == 0

	spec = obj.make_spec(overwrite=True)
	assert os.path.isfile(fn)
	assert os.stat(fn).st_size > 0


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


def test_SDSSObj_xid_datarelease():
	ra = 34.7489566190744
	dec = -4.04365821949229
	dir_obj = './testing/SDSSJ0218-0402/'


	run2ds = {12: 'v5_7_0', 13: 'v5_9_0'}

	for dr in [12, 13]:
		obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj, data_release=dr)
		assert obj.data_release == dr

		status = obj.load_xid(overwrite=True)

		assert status
		assert round(obj.xid['ra'], 4) == round(ra, 4)

		assert obj.run2d == run2ds[dr]


def test_SDSSObj_search_radius():

	search_radius = 5.*u.arcsec
	obj = sdssObj(ra=ra, dec=dec, dir_obj=dir_obj, search_radius=search_radius)

	assert obj.search_radius == search_radius


def test_SDSSObj_make_spec(obj_dirobj):
	obj = obj_dirobj
	status = obj.make_spec(overwrite=True)

	assert status
	assert os.path.isfile(obj.dir_obj + 'spec.fits')