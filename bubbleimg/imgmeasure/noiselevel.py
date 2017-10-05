# noiselevel.py
# ALS 2017/08/29

""" 
tool sets to measure the noise level of an image by fitting a Gaussian to the histogram of the image pixel values 
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.signal import medfilt

def getnoiselevel_gaussfit(data, fn_plot, toplot=True):
	"""
	return best fit gaussian noise level of the data, by fitting gaussian to
	the pixel value histogram. The histogram is medium filtered with a kernel size
	of 5 in order to mitigate the problem of spiky histogram of zero filled images. 

	Parameters:
	------
	data: 2d np array
		image
	toplot: bool

	Returns
	------
	sigma: float
	"""

	# set up
	nb = 1000 # number of bins

	histo, bin_edges = np.histogram(data.flatten(),bins=nb)
	std = np.std(data.flatten())

	bin_centers = bin_edges[0:-1]+np.diff(bin_edges)[0]/2.
	histo_med = medfilt(histo, kernel_size=5)

	p0 = [np.max(histo), 0., np.std(data.flatten())] # A, x0, sigma
	# p0 = [np.max(histo), std] # A, sigma
	coeff, var_matrix = curve_fit(gauss, bin_centers, histo_med, p0=p0)
	# sigma = np.abs(coeff[1])
	sigma = np.abs(coeff[2])

	if toplot:
		# plotting
		plt.clf()
		plt.axvline(0., color='grey', lw=0.5)
		plt.scatter(bin_centers, histo, label='data', color='#1f77b4', alpha=0.1)
		plt.scatter(bin_centers, histo_med, label='data med filtered', color='#1f77b4')
		plt.plot(bin_centers, gauss(bin_centers,*coeff), label='fit Gaussian sigma ='+'%.4e'%sigma, color='black')
		plt.legend()
		plt.xlabel('Intensity [img unit]')
		plt.ylabel('counts')
		plt.xlim(-2.5*std, 2.5*std)
		plt.savefig(fn_plot)
		plt.close()

	return sigma


def gauss(x, a, x0, sigma):
	return a*np.exp(-(x-x0)**2/(2*sigma**2))


def centered_gauss(x, a, sigma):
	x0 = 0.
	return a*np.exp(-(x-x0)**2/(2*sigma**2))

