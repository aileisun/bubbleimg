# test_extrap.py
# 06/07/2017 ALS 

"""
tool for integration
"""

import numpy as np
import matplotlib.pyplot as plt

from .. import extrap

def arr1_contains_arr2(arr1, arr2):

    return (max(arr1) >= max(arr2)) and (min(arr1) <= min(arr2))


def test_extrapolate():

	xs = np.logspace(3., 4., num=100)
	xs1 = xs[20:80]

	f = lambda x: 0.5*x + 1.

	ys1 = f(xs1)
	ys = f(xs)

	ys_test, xs_test = extrap.extrapolate(ys1, xs1, xs, polydeg=1, extbase_length=1000.)

	assert xs.size == xs_test.size
	assert np.all(np.around(xs, 3) == np.around(xs_test, 3))

	assert ys_test.size == ys.size
	assert np.all(np.around(ys, 3) == np.around(ys_test, 3))


def test_extrap_to_end_trailing():

	xs = np.logspace(3., 4., num=100)
	x_end = xs[-1]-1.
	xs1 = xs[:50]

	f = lambda x: 0.5*x + 1.

	ys1 = f(xs1)
	ys = f(xs)

	ys_test, xs_test = extrap.extrap_to_end(ys1, xs1, x_end, polydeg=1, extbase_length=1000.)

	assert xs.size == xs_test.size
	assert np.all(np.around(xs, 3) == np.around(xs_test, 3))

	assert ys_test.size == ys.size
	assert np.all(np.around(ys, 3) == np.around(ys_test, 3))


def test_extrap_to_end_heading():

	xs = np.logspace(3., 4., num=100)
	x_end = xs[0]+1.
	xs1 = xs[50:]

	f = lambda x: 0.5*x + 1.

	ys1 = f(xs1)
	ys = f(xs)
	ys_test, xs_test = extrap.extrap_to_end(ys1, xs1, x_end, polydeg=1, extbase_length=1000.)

	assert xs.size == xs_test.size
	assert np.all(np.around(xs, 3) == np.around(xs_test, 3))

	assert ys_test.size == ys.size
	assert np.all(np.around(ys, 3) == np.around(ys_test, 3))


def test_extrap_is_equalspace_log():
	xs1 = np.logspace(3., 4., num=100)

	assert extrap.is_equalspace_log(xs1)

	xs2 = np.linspace(3., 4., num=100)

	assert not extrap.is_equalspace_log(xs2)



def test_extrap_geomspace_stepsize():
	xs = np.logspace(3., 4., num=100)

	dlnx = extrap.get_dlnx(xs)

	xs_test = extrap.geomspace_stepsize(start=xs[0], stop=xs[-1], dlnx=dlnx)

	assert xs.size == xs_test.size
	assert np.all(xs == xs_test)

	# the new xs will go over stop a bit to cover it if num is not round number
	xs_test = extrap.geomspace_stepsize(start=xs[0], stop=xs[-1]-1., dlnx=dlnx)

	assert xs.size == xs_test.size
	assert np.all(xs == xs_test)


def test_extrap_linspace_stepsize():
	xs = np.linspace(1000., 10000., num=100)

	dx = extrap.get_dx(xs)

	xs_test = extrap.linspace_stepsize(start=xs[0], stop=xs[-1], dx=dx)

	assert xs.size == xs_test.size
	assert np.all(xs == xs_test)

	# the new xs will go over stop a bit to cover it if num is not round number
	xs_test = extrap.linspace_stepsize(start=xs[0], stop=xs[-1]-1., dx=dx)

	assert xs.size == xs_test.size
	assert np.all(xs == xs_test)


def test_extrap_extend_xs_to_end_trailing():

	# log
	xs = np.logspace(3., 4., num=100)

	xs_1h = xs[:50]
	x_end = xs[-1]-1.

	xs_test, __ = extrap.extend_xs_to_end(xs_1h, x_end)

	assert xs.size == xs_test.size
	assert np.all(np.around(xs, 3) == np.around(xs_test, 3))


def test_extrap_extend_xs_to_end_heading():

	# log
	xs = np.logspace(3., 4., num=100)

	xs_1h = xs[50:]
	x_end = xs[0]+1.

	xs_test, __ = extrap.extend_xs_to_end(xs_1h, x_end)

	assert xs.size == xs_test.size
	assert np.all(np.around(xs, 3) == np.around(xs_test, 3))