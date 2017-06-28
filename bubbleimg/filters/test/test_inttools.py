
import pytest
import numpy as np
import astropy.units as u

from .. import inttools
from ...obsobj import obsObj
from ... import spector

ra = 140.099341430207
dec = 0.580162492432517
dir_obj = '../../spector/test/testing/SDSSJ0920+0034/'


@pytest.fixture
def obj_dirobj():
	return obsObj(ra=ra, dec=dec, dir_obj = dir_obj)


@pytest.fixture
def spector1(obj_dirobj):
	obj = obj_dirobj
	return spector.Spector(obj=obj, survey_spec='boss', survey='hsc')


def test_inttools_calc_Fnu_in_band(spector1):
	s = spector1
	spec, ws = s.get_spec_ls()
	trans, ws_trans = s._get_norm_trans_func(band='i')

	Fnu = inttools.calc_Fnu_in_band_from_fl(fl=spec, ws=ws, trans=trans, ws_trans=ws_trans)


	assert (Fnu/u.Unit("erg s-1 cm-2 Hz-1")).unit == u.dimensionless_unscaled



def test_inttools_f():

	def f(x):
		return x

	# f is a function of x
	x = inttools.int_f_over_dlnx(f, x0=1., x1=2.)

	assert x == 1.


def test_inttools_arr():
	# same as above but integrate array instead of function

	def f(x):
		return x

	xs = np.linspace(1., 2., num=1001, endpoint=True)
	arr = f(xs)

	r0 = inttools.int_f_over_dlnx(f, x0=xs[0], x1=xs[-1])

	r1 = inttools.int_arr_over_dlnx(arr, xs)

	assert round(r0, 5) == round(r1, 5)



def test_inttools_arr_log():
	# same as above but integrate array instead of function

	def f(x):
		return x

	xs = np.logspace(0., 1., num=1001, endpoint=True)
	arr = f(xs)

	r0 = inttools.int_f_over_dlnx(f, x0=xs[0], x1=xs[-1])

	r1 = inttools.int_arr_over_dlnx(arr, xs)

	assert round(r0, 5) == round(r1, 5)



def test_inttools_f_and_arr():
	# int f1 * f2 dlnx
	def f1(x):
		return x

	def f2(x):
		return x + 1.

	def f(x):
		return f1(x)*f2(x)

	xs = np.logspace(0., 1., num=10001, endpoint=True)

	r0 = inttools.int_f_over_dlnx(f, x0=xs[0], x1=xs[-1])

	arr = f1(xs)
	r1 = inttools.int_arr_times_f_over_dlnx(arr, f2, xs)

	assert round(r0, 5) == round(r1, 5)


def test_inttools_arr_and_arr():
	# int f1 * f2 dlnx
	def f1(x):
		return x

	def f2(x):
		return x + 1.

	def f(x):
		return f1(x)*f2(x)

	xs1 = np.logspace(0., 1., num=10001, endpoint=True)
	xs2 = np.logspace(0., 1., num=8001, endpoint=True)

	r0 = inttools.int_f_over_dlnx(f, x0=xs1[0], x1=xs1[-1])

	arr1 = f1(xs1)
	arr2 = f2(xs2)
	r1 = inttools.int_arr_times_arr_over_dlnx(arr1, xs1, arr2, xs2)

	assert round(r0, 5) == round(r1, 5)


def test_inttools_arr1_contains_arr2():

	arr1 = np.linspace(0., 10.)
	arr2 = np.linspace(1., 9.)
	assert inttools.arr1_contains_arr2(arr1, arr2)

	arr1 = np.linspace(1., 10.)
	arr2 = np.linspace(1., 9.)
	assert inttools.arr1_contains_arr2(arr1, arr2)

	arr1 = np.linspace(2., 10.)
	arr2 = np.linspace(1., 9.)
	assert not inttools.arr1_contains_arr2(arr1, arr2)

	arr1 = np.linspace(0., 8.)
	arr2 = np.linspace(1., 9.)
	assert not inttools.arr1_contains_arr2(arr1, arr2)

