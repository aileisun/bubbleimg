# measure.py
# ALS 2016/04/11

"""
experimenting on various measurements that can be done on the img
"""

# from pylab import *
import os
import numpy as np
import matplotlib.pyplot as plt
# from astropy.io import fits

from astropy.table import Table,hstack,vstack
import astropy.units as u

import sys
sys.path.append( '/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/')
import measurenebula
reload(measurenebula)

import waveletdenoise

def main():

    """
    warning: the units of the measurements are not correct, should be in pixels, not arcsecs. 
    """
    # measure haar
    directory='haar/threshold_1sigma/'
    filename='img_waveletdenoised_soft_smooth.npy'

    run_measurements(directory,filename)

    # measure coif1
    directory='coif1/threshold_1sigma/'
    filename='img_waveletdenoised_soft.npy'

    run_measurements(directory,filename,iso_factor=1.)
    run_measurements(directory,filename,iso_factor=2.)

    
def run_measurements(directory,filename,iso_factor=1.):
    """
    make measurements on the file
    """
    # setup output
    fileout='measurements.csv'

    # read input
    img=np.load(directory+filename)
    params=Table.read(directory+'params.csv')

    sigma=params['sigma']
    isophotocut=sigma*iso_factor

    # make measurements
    tabiso, imgmos=measurenebula.measureISO(img,isophotocut=isophotocut,smoothing=1.,pixelsize=1.,useunit=False)
    params['isophotocut']=isophotocut
    tabout=hstack([params,tabiso])

    # write output
    if os.path.isfile(fileout):
        tabold=Table.read(fileout,format='csv')
        tabnew=vstack([tabold,tabout])
        tabnew.write(fileout,format='csv')
    else:
        tabout.write(fileout,format='csv')

    # plot corresponding clipped images that's being measured
    waveletdenoise.plot_img(img,'',vmin=sigma*iso_factor,vmax=30,tosave=False)
    plt.savefig(directory+filename[:-4]+'_clipped'+str(iso_factor)+'sigma.pdf')


def measure_moments(img,pixelsize=0.396,imgunit=u.Unit('erg s-1 cm-2 arcsec-2'),pixunit=u.Unit('arcsec')): 
    """
    PURPOSE: make measurements derived from image moments

    PARAMETERS: 
        img (np 2d array)         :  line intensity map (without unit)
        pixelsize=0.396           :  the image pixel size (without unit)
        imgunit=u.Unit('erg s-1 cm-2 arcsec-2')
                                  : image unit
        pixunit=u.Unit('arcsec')  : pixel size unit

    RETURN: 
        tabout: astropy table of 1 row with columns: 
             flux   : total integrated flux 
                      (default unit erg / (cm2 s)) 
             x_ctr  : x axis centroid relative to image center
                      (default unit arcsec)
             y_ctr  : y axis centroid relative to image center
                      (default unit arcsec) 
             sigma_major : length of the major axis (sigma)
                      (default unit arcsec)
             sigma_minor : length of the minor axis (sigma)
                      (default unit arcsec)
             theta  : orientation of the major axis, need to be tested
                      (unit degree)
             eccentricity : eccentricity of the ellipse 
                            sqrt(1-(minor/major)^2)

    DESCRIPTION: 
        0th - 2nd order image moments are calculated, and translated to 
        intuitave observables.  

        0th: 
            sum(I) :sum of pixel values   
        1st:
            sum(x * I): x centroids 
            sum(y * I): y centroids 
        2nd: 
            sum(x^2 * I), sum(y^2 * I), sum(xy * I): covariance matrix, 
            which can be translated to the size, eccentricity, and orientaiton of the ellipse. 

    REFERENCES: 
        https://en.wikipedia.org/wiki/Image_moment
        http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    """
    from skimage import measure
    from numpy import linalg
    from astropy.table import Table

    # sanity check, make sure units are consistent
    if hasattr(img, 'unit'):
        if img.unit != imgunit: raise ValueError('unit inconsistent')

    # calculate raw moments, m[ix,iy]
    m = measure.moments(img,order=2)
    # calculate centroids
    xc = m[1,0]/m[0,0]
    yc = m[0,1]/m[0,0]  
    # calculate central moments (mu in wiki), notice row = y, column=x.
    mc=measure.moments_central(img,cr=yc,cc=xc,order=2)
    # calculate normalized moments (mu' in wiki)
    mn=mc/mc[0,0]
    # make covariance matrix from 2nd order moments, in order to infer ellips
    cov=array([[mn[2,0],mn[1,1]],[mn[1,1],mn[0,2]]])
    # eigen decomposition of covariance matrix
    [lambda1,lambda2],vs=linalg.eig(cov)
    [lambda1,lambda2]=sort([lambda1,lambda2])[::-1]

    # infering measurements
    flux=m[0,0]*pixelsize**2 * (imgunit*pixunit**2)
    x_ctr = (xc-(img.shape[0]-1.)/2.)*pixelsize * (pixunit)
    y_ctr = (yc-(img.shape[0]-1.)/2.)*pixelsize * (pixunit)
    sigma_major, sigma_minor=sqrt([lambda1,lambda2])*pixelsize * (pixunit)
    theta=-degrees(arctan(vs[1,0]/vs[0,0]))*u.Unit('degree')
    eccentricity=sqrt(1-lambda2/lambda1)
    
    # prep output
    tabout=Table()
    for var in ['flux','x_ctr','y_ctr','sigma_major','sigma_minor','theta','eccentricity']: 
        if isinstance(eval(var),astropy.units.quantity.Quantity):
            tabout[var]=[eval(var).value]
            tabout[var].unit=eval(var).unit
        else: 
            tabout[var]=[eval(var)]

    return tabout




def overplot_moment_ellips(img,tab_moments,pixelsize=0.396):
    """
    PURPOSE: 
        plot the corresponding ellipse on top of the image given moment measurements
    INPUT: 
        img: image array
        tab_moments: table containing momoent measurements with columns including: 'x_ctr','y_ctr','sigma_major','sigma_minor','theta'
        pixelsize=0.396

    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse    

    # unit conversion from arcsec to pix
    for col in ['x_ctr','y_ctr','sigma_major','sigma_minor']: # sanity check
        if tab_moments[col].unit != u.Unit('arcsec'): 
            raise NameError('check conversion')
    xc = tab_moments['x_ctr'][0]/pixelsize+(img.shape[0]-1.)/2.
    yc = tab_moments['y_ctr'][0]/pixelsize+(img.shape[0]-1.)/2.
    smajor = tab_moments['sigma_major'][0]/pixelsize
    sminor = tab_moments['sigma_minor'][0]/pixelsize
    theta = tab_moments['theta'][0]

    # define ellipse
    kwrg = {'facecolor':'none', 'edgecolor':'black', 'alpha':1, 'linewidth':2}
    ellip = Ellipse(xy=[xc,yc], width=2*smajor, height=2*sminor, angle=theta, **kwrg)
    
    # plotting image and ellipse
    plt.clf()
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    ax.imshow(img,interpolation='nearest',origin='lower left',aspect='equal',cmap='viridis',extent=[0,img.shape[0],0,img.shape[0]])
    ax.plot(xc,yc,marker='x',ms=5,mew=2,color='black')
    ax.add_artist(ellip)
    ax.set_xlim(0,img.shape[0])
    ax.set_ylim(0,img.shape[0])


