# getFilterResponseFunc.py
# ALS 2016/05/02

import numpy as np
from astropy.io import fits
import astropy.units as u

from astropy.table import Table

filterwavelengths={'u': 3551.*u.AA, 'g': 4686.*u.AA, 'r': 6166.*u.AA, 'i': 7480.*u.AA, 'z': 8932.*u.AA}


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

