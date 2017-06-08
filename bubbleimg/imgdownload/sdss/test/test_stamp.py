# test_stamp.py
# ALS 2017/03/17

"""
to be used with pytest

test sets for stamp
"""

import numpy as np
from .. import stamp


def test_cutstampImage_upperright():
	"""
	test that cutstampImage correctly cut stamp and fill in blanc values
	"""

	image = np.ones([4, 4])
	xcenter = 3
	ycenter = 3
	xwidth = 3
	ywidth = 3
	stamp_truth = np.array([[1, 1, 0,], [1, 1, 0,], [0, 0, 0,], ])

	stmp = stamp.cutstampImage(image, xcenter, ycenter, xwidth, ywidth)

	assert np.all(stmp == stamp_truth)


def test_cutstampImage_lowerleft():
	"""
	test that cutstampImage correctly cut stamp and fill in blanc values
	"""

	image = np.ones([4, 4])
	xcenter = 0
	ycenter = 0
	xwidth = 3
	ywidth = 3
	stamp_truth = np.array([[0, 0, 0,], [0, 1, 1,], [0, 1, 1,], ])

	stmp = stamp.cutstampImage(image, xcenter, ycenter, xwidth, ywidth)

	assert np.all(stmp == stamp_truth)

