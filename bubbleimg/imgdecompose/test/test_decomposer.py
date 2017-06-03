
import pytest
import shutil
import os

from ...obsobj import obsObj

from .. import decomposer

ra = 140.099341430207
dec = 0.580162492432517
dir_parent = './testing/'
dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'test_verification_data/SDSSJ0920+0034/'
bandline = 'i'
bandconti = 'r'
survey = 'hsc'

@pytest.fixture(scope="module", autouse=True)
def setUp_tearDown():
	""" rm ./testing/ and ./test2/ before and after testing"""

	# setup
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)

	os.makedirs(dir_parent)
	shutil.copytree(dir_verif, dir_obj)

	yield
	# tear down
	if os.path.isdir(dir_parent):
		shutil.rmtree(dir_parent)


@pytest.fixture
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


def test_decomping(obj_dirobj):
	obj = obj_dirobj

	d = decomposer.Decomposer(obj=obj, survey='hsc')
	assert d.survey == survey

	d.make_spec_SEDs(overwrite=False)

	assert os.path.isfile(obj.dir_obj+'sed.csv')
	assert os.path.isfile(obj.dir_obj+'sed_conti.csv')
	assert os.path.isfile(obj.dir_obj+'sed_line.csv')


	d.make_contisub(bandline=bandline, bandconti=bandconti, z=None, overwrite=False)

	assert os.path.isfile(obj.dir_obj+'stamp-{0}-{1}conti.fits'.format(bandline, bandconti))

	d.make_linemap(line='oiii5007', bandline=bandline, bandconti=bandconti, z=None, overwrite=False)

	assert os.path.isfile(obj.dir_obj+'stamp-lOIII5008_I.fits')



