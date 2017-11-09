# test_getzrange_batch.py
# ALS 2017/06/28

import numpy as np

import pytest
from .. import getzrange_batch
import imp
imp.reload(getzrange_batch)


def test_determine_non_excluded_zrange():

	in_zmin = 0.
	in_zmax = 10.

	# case exclude lower end
	ex_zmin = -1.
	ex_zmax = 2.

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax)

	assert z0 == ex_zmax
	assert z1 == in_zmax


	# case exclude higher end
	ex_zmin = 9
	ex_zmax = 11.

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax)

	assert z0 == in_zmin
	assert z1 == ex_zmin


	# case exclude middle part
	ex_zmin = 1
	ex_zmax = 9

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax)


	assert np.all(z0 == [in_zmin, ex_zmax])
	assert np.all(z1 == [ex_zmin, in_zmax])


	# case excluded entirely
	ex_zmin = -1
	ex_zmax = 11

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax)

	assert np.isnan(z0)
	assert np.isnan(z1)



	# case no exclusion
	ex_zmin = -5
	ex_zmax = -1

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax)

	assert z0 == in_zmin
	assert z1 == in_zmax



def test_determine_non_excluded_zrange_wrong_order():
	# order should not matter
	in_zmin = 0.
	in_zmax = 10.

	# case exclude lower end
	ex_zmin = -1.
	ex_zmax = 2.

	z0, z1 = getzrange_batch.determine_non_excluded_zrange(in_zmin=in_zmax, in_zmax=in_zmin, ex_zmin=ex_zmax, ex_zmax=ex_zmin)

	assert z0 == ex_zmax
	assert z1 == in_zmax
