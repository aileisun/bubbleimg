# filtertools.py
# ALS 2016/05/02

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u

from astropy.table import Table

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy.optimize import fsolve


filterwavelengths={'u': 3551.*u.AA, 'g': 4686.*u.AA, 'r': 6166.*u.AA, 'i': 7480.*u.AA, 'z': 8932.*u.AA}


def findFilterBounday(threshold=0.6,toplot=True):
    """
    PURPOSE: write file filterboundary.txt to record filter boundary defined by 80% of maximum throughput
    PARAMETERS:
            threshold=0.6 (float)
            toplot=True
    """

    # fileout
    fileout='filterboundary_'+'%.1f'%threshold+'.txt'

    tabout=Table([[],[],[],],names=('band','w1','w2'),dtype=('string','int','int'))

    for band in ['u','g','r','i','z']:
        # readin filter function
        spec,lcoord=getFilterResponseFunc(band=band)


        spl=UnivariateSpline(lcoord, spec/max(spec)-threshold,s=0)
        roots=spl.roots()
        # print roots

        # solve for root
        f = interp1d(lcoord, spec/max(spec),kind='linear',bounds_error=False,fill_value=0.)
        tosolve=lambda x: f(x)-threshold

        x1=fsolve(tosolve,x0=roots[0])
        x2=fsolve(tosolve,x0=roots[-1])

        tabout.add_row([band,x1,x2])

        if toplot:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(lcoord,spec/max(spec),color='black')
            ax.axhline(threshold)
            ax.axvline(x1)
            ax.axvline(x2)
            ax.set_xlim(min(lcoord),max(lcoord))
            # raw_input("Enter...")
            fig.savefig(band+'_'+'%.1f'%threshold+'.pdf')

    tabout.write(fileout,format='ascii.fixed_width',delimiter='')


def intRoverldl(band='r'):
    """
    For a specified band, calcualte and return the integral 
        s=int{R(l)/l dl}, 
    where R is the filter transmission funciton. 
    s is a dimensionless quantity. 
    """
    R,l=getFilterResponseFunc(band=band)
    # sanity check - dl constant
    if len(np.unique(np.diff(l)))==1: dl=np.diff(l)[0]
    else: raise ValueError("Filter response function not equally spaced in lambda")
    # integral
    s=0.
    for i in range(len(R)):
        s=s+R[i]/l[i]*dl
    return s



def getlocalpath():
    """
    return path to filter direcotry
    """
    import os
    import sys 
    path=os.path.dirname(sys.modules[__name__].__file__)
    if path == '': path ='.'
    return path+'/'



def accessFile(filename='OIIIredshiftrange0.6.txt'):
    """
    access files in filters/ such as:
        HaNIIredshiftrange0.2.txt
        OIIIredshiftrange0.6.txt
    """
    pathlocal=getlocalpath()
    return Table.read(pathlocal+filename,format='ascii')
    


def getFilterResponseFunc(band='r'):
    """
    PURPOSE: 
        Returning filter response funciton given band, which includes extinction through 
        an airmass of 1.3 at Apache Point Observatory. Note that these are not complete
        filter curves, as they do not include the full system response from atmosphere 
        to detector.

        See http://www.sdss.org/instruments/camera/#Filters
    PAREMATER: 
        band='r' (string) band has to be one of ['u','g','r','i','z']

    RETURN: 
        respt      (array) response funciton
        wavelength (array) [AA]
    """
    # read filter funciton
    filename=getlocalpath()+'filter_curves.fits'
    hdulist=fits.open(filename)

    # selecting filter funciton
    bands=np.array(['u','g','r','i','z'])
    ib=np.where(bands==band)[0]

    if ib.size != 1: raise ValueError("input band not recognized")
    else: ib=ib[0]
    if not hdulist[ib+1].header['EXTNAME']==band.capitalize(): raise ValueError("Error in matching band")

    R=hdulist[ib+1].data['respt']
    l=hdulist[ib+1].data['wavelength']
    return R, l

