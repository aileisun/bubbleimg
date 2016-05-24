# noiselevel.py
# ALS 2016/04/11


"""
measuring noise levels
"""

from pylab import *

from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.table import Table


# def main(directory='./'):
#     """
#     estimate, write, and plot sigma
#     """
#     # setting
#     scaling = 1.e15

#     # read in data
#     filename=directory+'stamp-lOIII5008_I.fits'
#     fileout=directory+'noiselevel.csv'
#     hdu=fits.open(filename)
#     img=hdu[0].data*scaling # scale up the data

#     sigma=getnoiselevel(img,toplot=True)
#     peaksnr=np.max(img)/sigma
#     print 'sigma, snr',sigma,peaksnr

#     tsummary=Table([[filename],[sigma],[peaksnr]],names=['filename','sigma',['peaksnr']])
#     tsummary.write(fileout,format='csv')


def getnoiselevel(data,ranges=(-30,30),toplot=False):
    """
    return best fit gaussian error given data
    """
    nb=100

    histo,bin_edges=np.histogram(data.flatten(),bins=nb,range=ranges)

    bin_centers = bin_edges[0:-1]+np.diff(bin_edges)[0]/2.

    p0 = [4000., 0, 1.] # A, x0, sigma
    coeff, var_matrix = curve_fit(gauss, bin_centers, histo, p0=p0)
    coeff[1:3]=coeff[1:3]
    sigma=coeff[2]

    if toplot:
        # plotting
        plt.clf()
        plt.plot(bin_centers,histo,label='data')
        plt.plot(bin_centers,gauss(bin_centers,*coeff),label='fit Gaussian sigma ='+'%.4e'%sigma)
        plt.legend()
        plt.xlabel('Intensity [img unit]')
        plt.ylabel('counts')
        plt.savefig('noiselevel.pdf')

    return sigma

def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


