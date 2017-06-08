# filtertools.py
# ALS 2016/05/02

import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u

import astropy.table as at

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

import glob

from surveysetup import surveybands
from surveysetup import waverange

def getFilterCentroids(band='u', survey='sdss', withunit=True):
    """
    Return the table of filter centroids given the survey
    """
    fc = accessFile(filename='filtercentroid.txt', survey=survey, joinsurveys=True)

    if withunit: 
        unit = u.AA
    else:
        unit = 1.

    return fc['w'][fc['band']==band][0]*unit


def writeFilterCentroids(survey='sdss'):
    """
    Purpose: write file filtercentroid.txt to record filter centroid defined as transmission weighted average wavelength

    Params
    ----------
    survey = 'sdss'

    Return
    ---------
    None

    Output
    ----------
    filtercentroid.txt
    """
    localpath = getlocalpath()
    fileout
    bands = surveybands[survey]

    tabout=at.Table([[],[],],names=('band','w'),dtype=('string','int'))

    for band in bands:
        # readin filter function
        spec,lcoord = getFilterResponseFunc(band=band, survey=survey)
        w = np.average(lcoord, weights=spec)

        tabout.add_row([band, w])

    tabout.write(fileout,format='ascii.fixed_width',delimiter='')


def writeFilterBoundaries(threshold=0.6, toplot=True, survey='sdss'):
    """
    PURPOSE: write file filterboundary.txt to record filter boundary defined by threshold*100% of maximum throughput

    PARAMETERS:
            threshold=0.6: (float)
            toplot=True: (bool)
            survey='sdss': (string)
                    one of the following: sdss, hsc, ukirt
    Return
    ---------
    None

    Output
    ----------
    filterboundary.txt
    """

    # fileout
    localpath = getlocalpath()
    fileout
    bands = surveybands[survey]


    tabout=at.Table([[],[],[],],names=('band','w1','w2'),dtype=('string','int','int'))

    for band in bands:
        # readin filter function
        spec,lcoord=getFilterResponseFunc(band=band, survey=survey)

        x = lcoord
        y = spec/max(spec)-threshold

        spl=UnivariateSpline(x, y,s=0)
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



def getFilterResponseFunc(band='r', survey='sdss'):
    """
    PURPOSE: 
        Returning filter response funciton given band, which includes extinction through 
        an airmass of 1.3 at Apache Point Observatory. Note that these are not complete
        filter curves, as they do not include the full system response from atmosphere 
        to detector.

        See http://www.sdss.org/instruments/camera/#Filters
    PAREMATER: 
        band='r' (string) band has to be one of ['u','g','r','i','z','y','j','h','k'] depending on the survey
        survey = 'sdss' (string): among sdss, hsc, ukirt

    RETURN: 
        resp      (array) response funciton
        wavelength (array) [AA]
    """
    # read filter funciton

    if survey == 'sdss': 
        filename=getlocalpath()+survey+'/'+'filter_curves.fits'
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
    elif survey == 'hsc': 
        filename=getlocalpath()+survey+'/'+'filter_curves/'+'HSC-'+band+'.txt'
        tab = at.Table.read(filename, format='ascii')
        R = np.array(tab['col2'])
        l = np.array(tab['col1']*10)

        if l[1] < l[0]: 
            R = R[::-1]
            l = l[::-1]
        return R, l
    elif survey == 'ukirt':
        filename=getlocalpath()+survey+'/'+'filter_curves/'+'ukirt-'+band+'.txt'
        tab = at.Table.read(filename, format='ascii', comment='#')
        R = np.array(tab['R'])
        l = np.array(tab['l']*10)
        if l[1] < l[0]: 
            R = R[::-1]
            l = l[::-1]
        return R, l
    elif survey == 'cfht':
        filename=getlocalpath()+survey+'/'+'filter_curves/'+'cfht-'+band+'.txt'
        tab = at.Table.read(filename, format='ascii', comment='#')
        l = np.array(tab['col1'])
        R = np.array(tab['col2'])
        if l[1] < l[0]: 
            R = R[::-1]
            l = l[::-1]
        return R, l
    else: 
        raise NameError('survey name not recognized. ')


def intRoverldl(band='r', survey='sdss'): 
    """
    For a specified band, calcualte and return the integral 
        s=int{R(l)/l dl}, 
    where R is the filter transmission funciton. 
    s is a dimensionless quantity. 
    """
    R, l = getFilterResponseFunc(band=band, survey=survey)

    # sanity check - dl constant
    if len(np.unique(np.diff(l)))==1: 
        dl=np.diff(l)[0]
    else: 
        raise ValueError("Filter response function not equally spaced in lambda")
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



def accessFile(filename='zranges_band_wOIII_nHaNIISII.txt', survey='sdss', joinsurveys=True):
    """
    access files in filters/ such as:
        HaNIIredshiftrange0.2.txt
        OIIIredshiftrange0.6.txt

    PARAMS
    ---------
    filename='OIIIredshiftrange0.6.txt'
    survey='sdss': (string)
            for joined surveys use '-', e.g., hsc-ukirt
    joinsurveys=True (string)
            if true then if survey/filename does not exist, for joined surveys go into each directories/files and merge file
    """

    filepath = getlocalpath()+survey+'/'+filename
    if os.path.isfile(filepath):
        tab = at.Table.read(filepath, format='ascii')
    elif ('-' in survey) and joinsurveys: 
        surveys = survey.split('-')
        tab = at.Table()
        for s in surveys:
            filepath = getlocalpath()+s+'/'+filename
            tabnew = at.Table.read(filepath, format='ascii')
            tab = at.vstack([tab, tabnew])
    else:
        raise NameError('File does not exist')

    return tab

    
def accessTabZranges(lineconfig = 'wOIII_nHaNIISII', survey='sdss', joinsurveys=True):
    """
    access zranges table

    PARAMS
    ----------
    lineconfig = 'wOIII_nHaNIISII'
        if 'all' than all the tables for all the line cofigs available in the survey is combined and returned
    surver = 'sdss'
    joinsurveys=True
    """

    if lineconfig !='all':
        filenames = ['zranges_band_'+lineconfig+'.txt']
        tab = accessFile(filename=f, survey=survey, joinsurveys=joinsurveys)

    else:
        dir_survey = getlocalpath()+survey+'/'
        filenames = glob.glob(dir_survey+'zranges_band_*.txt')
        tab = at.Table()
        for f in filenames:
            tabnew = at.Table.read(f, format='ascii')
            tab = at.vstack([tab, tabnew])

        print filenames

    return tab