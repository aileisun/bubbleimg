
from .. import filtertools

from ...spector import inttools
from ..surveysetup import surveybands

# import pytest

# dir_test = 'testing/'

# @pytest.fixture(scope="module", autouse=True)
# def setUp_tearDown():
# 	""" rm ./testing/ and ./test2/ before and after testing"""

# 	# setup
# 	if os.path.isdir(dir_test):
# 		shutil.rmtree(dir_test)
# 	os.makedirs(dir_test)

# 	yield
# 	# tear down
# 	if os.path.isdir(dir_test):
# 		shutil.rmtree(dir_test)


def test_filtertools_calcintresdlnl():

	band = 'r'
	survey = 'hsc'

	a = filtertools.calc_int_response_dlnl(band=band, survey=survey)

	trans, ws_trans = filtertools.getFilterResponseFunc(band=band, survey=survey)

	a_test = inttools.int_arr_over_dlnx(arr=trans, xs=ws_trans)

	assert a == a_test


def test_filtertools_getintresdlnl():

	band = 'r'
	survey = 'hsc'

	a = filtertools.get_int_response_dlnl(band=band, survey=survey)

	trans, ws_trans = filtertools.getFilterResponseFunc(band=band, survey=survey)

	a_test = inttools.int_arr_over_dlnx(arr=trans, xs=ws_trans)

	assert round(a, 3) == round(a_test, 3)


def test_filtertools_calcNormTransFunc():

	for survey in ['sdss', 'hsc']:
		for band in surveybands[survey]:

			trans, ws = filtertools.calcNormTransFunc(band=band, survey=survey)

			inttransdlnl = inttools.int_arr_over_dlnx(arr=trans, xs=ws)
			assert  round(inttransdlnl, 5) == 1.


def test_filtertools_getNormTransFunc():

	for survey in ['sdss', 'hsc']:
		for band in surveybands[survey]:

			trans, ws = filtertools.getNormTransFunc(band=band, survey=survey)

			inttransdlnl = inttools.int_arr_over_dlnx(arr=trans, xs=ws)
			assert  round(inttransdlnl, 5) == 1.


def test_filtertools_getNormTrans():

	f1 = filtertools.getNormTrans(l=5520.0, band='r', survey='hsc')
	assert f1 == 3.4454694768050853

	f2 = filtertools.getNormTrans(l=5530.0, band='r', survey='hsc')
	f3 = filtertools.getNormTrans(l=5540.0, band='r', survey='hsc')

	assert f2 == (f1+f3)/2.

	for survey in ['sdss', 'hsc']:
		for band in surveybands[survey]:

			trans, ws = filtertools.getNormTransFunc(band=band, survey=survey)

			inttransdlnl = inttools.int_arr_over_dlnx(arr=trans, xs=ws)
			assert  round(inttransdlnl, 5) == 1.

