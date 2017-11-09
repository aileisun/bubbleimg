
import numpy as np
import matplotlib.pyplot as plt


from . import polytools
from . import plottools

def main():
	img = np.array([[0, 0, 0,], [0, 1, 0], [0, 0, 0]])
	threshold = 0.5

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	area = polytools.NetPolygonsArea(contours)

	fig, ax = plottools.plot_img(img, colorlabel='')
	plottools.overplot_contours(ax, contours)



	img = np.array([[0, 0, 0, 0,], [0, 1, 0, 0,], [0, 1, 0, 0,], [0, 0, 0, 0,], ])
	threshold = 0.9

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	area = polytools.NetPolygonsArea(contours)
	print(area)

	fig, ax = plottools.plot_img(img, colorlabel='')
	plottools.overplot_contours(ax, contours)




	img = np.array([[0, 0, 0, 0, 0,], [0, 1, 1, 1, 0, ], [0, 1, 0, 1, 0, ], [0, 1, 1, 1, 0, ], [0, 0, 0, 0, 0,], ])
	threshold = 0.5

	contours = polytools.find_contours(img, threshold, tocomplete=True)

	area = polytools.NetPolygonsArea(contours)
	print(area)


	fig, ax = plottools.plot_img(img, colorlabel='')
	plottools.overplot_contours(ax, contours)
