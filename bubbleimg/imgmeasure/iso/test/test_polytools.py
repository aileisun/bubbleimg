import pytest

import numpy as np
import matplotlib.pyplot as plt

from .. import polytools
from .. import plottools

# this is just a simple test of the basic tests in polytools. This is by no mean an extensive test set. 

dir_test = 'testing/'

img_dot = np.array([[0, 0, 0,], [0, 1, 0], [0, 0, 0]])
img_donut = np.array([[0, 0, 0, 0, 0,], [0, 1, 1, 1, 0, ], [0, 1, 0, 1, 0, ], [0, 1, 1, 1, 0, ], [0, 0, 0, 0, 0,], ])
img_eight = np.array([[0, 0, 0, 0, 0,], [0, 1, 1, 1, 0, ], [0, 1, 0, 1, 0, ], [0, 1, 1, 1, 0, ], [0, 1, 0, 1, 0, ], [0, 1, 1, 1, 0, ], [0, 0, 0, 0, 0,], ])

thrshld = 0.5


def test_find_contours_one():

	img = img_dot
	threshold = thrshld

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	assert len(contours) == 1.

	print(("real area: ", np.sum(img), "cntr area: ", polytools.NetPolygonsArea(contours)))

	assert about_the_same_area(contours, np.sum(img), uncertainty=0.)


def test_find_contours_donut():
	img = img_donut
	threshold = thrshld

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	assert len(contours) == 2

	highcontour = polytools.select_highcontours(contours)
	assert len(highcontour) == 1

	lowcontour = polytools.select_lowcontours(contours)
	assert len(lowcontour) == 1

	assert polytools.isinside_polygon(lowcontour[0], highcontour[0])

	print(("real area: ", np.sum(img), "cntr area: ", polytools.NetPolygonsArea(contours)))

	assert about_the_same_area(contours, np.sum(img), uncertainty=0.)


def test_find_contours_eight():
	img = img_eight
	threshold = thrshld

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	assert len(contours) == 3

	highcontour = polytools.select_highcontours(contours)
	assert len(highcontour) == 1

	lowcontour = polytools.select_lowcontours(contours)
	assert len(lowcontour) == 2

	assert polytools.isinside_polygon(lowcontour[0], highcontour[0])

	print(("real area: ", np.sum(img), "cntr area: ", polytools.NetPolygonsArea(contours)))

	assert about_the_same_area(contours, np.sum(img), uncertainty=0.)


def test_find_centercontours():

	# one
	img = img_dot
	threshold = thrshld

	xc = 1.
	yc = 1.

	contours = polytools.find_centercontours(img, threshold, xc, yc, radius=0.)
	assert len(contours) == 1

	contours = polytools.find_centercontours(img, threshold, xc, yc, radius=1.)
	assert len(contours) == 1

	xc = 0.
	yc = 0.

	contours = polytools.find_centercontours(img, threshold, xc, yc, radius=0.)
	print(contours)
	assert len(contours) == 0

	contours = polytools.find_centercontours(img, threshold, xc, yc, radius=2.)
	assert len(contours) == 1.

	# donut
	img = img_donut
	xc = 2.
	yc = 2.

	contours = polytools.find_centercontours(img, threshold, xc, yc, radius=0.)
	highcontours = polytools.select_highcontours(contours)
	assert len(highcontours) == 1

	lowcontours = polytools.select_lowcontours(contours)
	assert len(lowcontours) == 1


def test_find_largecontours():

	# one
	threshold = thrshld

	# =======================
	areallimit = 1. 

	contours = polytools.find_largecontours(img_dot, threshold, areallimit)
	assert len(contours) == 1.

	contours = polytools.find_largecontours(img_donut, threshold, areallimit)
	assert len(contours) == 2.

	contours = polytools.find_largecontours(img_eight, threshold, areallimit)
	assert len(contours) == 3.


	# =======================
	areallimit = 8. 

	contours = polytools.find_largecontours(img_dot, threshold, areallimit)
	assert len(contours) == 0.

	contours = polytools.find_largecontours(img_donut, threshold, areallimit)
	assert len(contours) == 2.

	contours = polytools.find_largecontours(img_eight, threshold, areallimit)
	assert len(contours) == 3.


	# =======================
	areallimit = 9. 

	contours = polytools.find_largecontours(img_dot, threshold, areallimit)
	assert len(contours) == 0.

	contours = polytools.find_largecontours(img_donut, threshold, areallimit)
	assert len(contours) == 0.

	contours = polytools.find_largecontours(img_eight, threshold, areallimit)
	assert len(contours) == 3.


def tests_ShapeParamsTab_from_contours():

	threshold = thrshld
	xc, yc = 2, 2
	contours = polytools.find_contours(img_donut, threshold)

	tab = polytools.ShapeParamsTab_from_contours(contours, xc, yc)

	print("")
	print(tab)
	assert len(tab) == 1

	assert about_the_same_area(contours, tab['area_pix'], uncertainty = 0.) 

	assert np.absolute(tab['dmax_pix'] - 4.24) < 0.8

	assert tab['dmax_pix'] == 2.*tab['rmax_pix']


def test_plot():
	img = img_donut
	threshold = thrshld

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	fig, ax = plottools.plot_img(img, colorlabel='')
	plottools.overplot_contours(ax, contours)
	fig.savefig(dir_test+'img_donut.pdf')


def about_the_same_area(contours, area, uncertainty = 0.):
	""" check if the conoturs enclsoe about the same area as area """

	carea = polytools.NetPolygonsArea(contours)

	return np.absolute(carea - area) <= uncertainty
