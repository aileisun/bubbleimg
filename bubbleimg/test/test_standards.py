from .. import standards
import numpy as np


def test_get_img_xycenter():

	img = np.ndarray([3, 3])
	xc, yc = standards.get_img_xycenter(img, center_mode='n/2')

	assert xc == 1
	assert yc == 1

	img = np.ndarray([4, 4])
	xc, yc = standards.get_img_xycenter(img, center_mode='n/2')

	assert xc == 2
	assert yc == 2

	img = np.ndarray([4, 4])
	xc, yc = standards.get_img_xycenter(img, center_mode='n/2-1')

	assert xc == 1
	assert yc == 1


	img = np.ndarray([3, 4])
	xc, yc = standards.get_img_xycenter(img, center_mode='n/2')

	assert xc == 1
	assert yc == 2