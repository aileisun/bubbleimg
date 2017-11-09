# filtertools.py
# ALS 2016/05/02

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u

import astropy.table as at

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

import glob

from .surveysetup import surveybands
from .surveysetup import waverange

from . import inttools


def getlocalpath():
    """
    return path to filter direcotry
    """
    path = os.path.dirname(sys.modules[__name__].__file__)
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
        raise NameError('File does not exist: '+filepath)

    return tab


def isFile(filename='zranges_band_wOIII_nHaNIISII.txt', survey='sdss'):
    """
    return whether there is such a file under survey 
    """
    filepath = getlocalpath()+survey+'/'+filename

    return os.path.isfile(filepath)

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

        R = hdulist[ib+1].data['respt']
        l = hdulist[ib+1].data['wavelength']
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
    fileout = localpath+survey+'/filtercentroid.txt'
    bands = surveybands[survey]

    tabout = at.Table([[], [], ], names=('band', 'w'), dtype=('string', 'int'))

    for band in bands:
        # readin filter function
        spec, ws = getFilterResponseFunc(band=band, survey=survey)
        w = np.average(ws, weights=spec)

        tabout.add_row([band, w])

    tabout.write(fileout, format='ascii.fixed_width', delimiter='')



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
    print("[filtertools] running writeFilterBoundaries()")
    localpath = getlocalpath()
    if np.log10(threshold) > -1:
        fileout = localpath+survey+'/'+get_fn_FilterBoundaries(threshold)
    else:
        fileout = localpath+survey+'/'+get_fn_FilterBoundaries(threshold)

    bands = surveybands[survey]


    tabout = at.Table([[], [], [], ], names=('band', 'w1', 'w2'), dtype=('string', 'int', 'int'))

    for band in bands:
        # readin filter function
        spec, ws = getFilterResponseFunc(band=band, survey=survey)

        x = ws
        y = spec/max(spec) - threshold

        spl = UnivariateSpline(x, y, s=0)
        roots = spl.roots()

        # solve for root
        f = interp1d(ws, spec/max(spec), kind='linear', bounds_error=False, fill_value=0.)
        tosolve = lambda x: f(x)-threshold

        x1 = fsolve(tosolve, x0=roots[0])
        x2 = fsolve(tosolve, x0=roots[-1])

        tabout.add_row([band, x1, x2])

        if toplot:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(ws, spec/max(spec), color='black')
            ax.axhline(threshold)
            ax.axvline(x1)
            ax.axvline(x2)
            ax.set_xlim(min(ws) ,max(ws))
            # raw_input("Enter...")
            fig.savefig(band+'_'+'%.1f'%threshold+'.pdf')

    tabout.write(fileout, format='ascii.fixed_width', delimiter='')


def get_fn_FilterBoundaries(threshold):
    """ return the file name of filterbounday given threshold """
    if np.log10(threshold) > -1:
        return 'filterboundary_'+'%.1f'%threshold+'.txt'
    else:
        return 'filterboundary_'+'%.2f'%threshold+'.txt'


def getFilterBoundaries(threshold=0.6, band='u', survey='sdss', withunit=True):
    """
    Return the filter boundary of a band in a survey given threshold
    """
    fn = get_fn_FilterBoundaries(threshold)

    if not isFile(filename=fn, survey=survey):
        writeFilterBoundaries(threshold=threshold, toplot=False, survey=survey)

    tab = accessFile(filename=fn, survey=survey, joinsurveys=False)

    i = [tab['band']==band]
    w1 = tab['w1'][i][0]
    w2 = tab['w2'][i][0]

    if withunit: 
        unit = u.AA
    else:
        unit = 1.

    return w1*unit, w2*unit


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

        print(filenames)

    return tab


def calc_int_response_dlnl(band='r', survey='sdss'):
    """ calculate int{ T(l)  d lnl } """

    trans, ws_trans = getFilterResponseFunc(band=band, survey=survey)

    a = inttools.int_arr_over_dlnx(arr=trans, xs=ws_trans)

    return a


def write_int_response_dlnl(survey='sdss'):
    """ write file  intresdlnl.txt  for the survey that records the values of int res dlnl for each band """
    localpath = getlocalpath()
    fileout = localpath+survey+'/intresdlnl.txt'

    bands = surveybands[survey]

    tabout = at.Table([[], [], ], names=('band', 'inttransdlnl'), dtype=('string', 'float'))

    for band in bands:
        # readin filter function

        trans, ws_trans = getFilterResponseFunc(band=band, survey=survey)

        result = inttools.int_arr_over_dlnx(arr=trans, xs=ws_trans)

        tabout.add_row([band, result])

    tabout.write(fileout, format='ascii.fixed_width', delimiter='', overwrite=True)


def get_int_response_dlnl(band='r', survey='sdss'):
    """ return precalcualted int trans dlnl of the specified band and survey """
    localpath = getlocalpath()
    fn = localpath+survey+'/intresdlnl.txt'

    if not os.path.isfile(fn):
        write_int_response_dlnl(survey=survey)

    tab = accessFile(filename='intresdlnl.txt', survey=survey, joinsurveys=False)

    return tab['inttransdlnl'][tab['band']==band][0]


def calcNormTransFunc(band='r', survey='sdss'):
    """ 
    calculate normalized filter transmission function from raw response funciton
    the transmission function is normalized such that int {T{l} dlnl} = 1.

    Return
    ------
    res (array)
    ws  (array)
    """

    res, ws = getFilterResponseFunc(band=band, survey=survey)
    intresdlnl = get_int_response_dlnl(band=band, survey=survey)

    return res/intresdlnl, ws


def writeNormTransFunc(survey='sdss'):
    """ 
    write directory under the survey directory

    normtrans/
        band1.csv
        band2.csv
        ...
    that records the normalized filter transmission function that has int {T{l} dlnl} = 1.
    """
    localpath = getlocalpath()
    dir_out = localpath+survey+'/normtrans/'

    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)

    for band in surveybands[survey]:
        fp = dir_out + band+'.csv'

        trans, ws = calcNormTransFunc(band=band, survey=survey)

        tabout = at.Table([ws, trans], names=('ws', 'trans'), dtype=('float', 'float'))

        tabout.write(fp, format='ascii.csv', overwrite=True)


def getNormTransFunc(band='r', survey='sdss'):
    """ return the normalized filter transmission function that is precalculated and stored in normtrans/ """
    localpath = getlocalpath()
    fn = localpath+survey+'/normtrans/'+band+'.csv'

    if not os.path.isfile(fn):
        writeNormTransFunc(survey=survey)

    tab = at.Table.read(fn, format='ascii.csv')

    ws = np.array(tab['ws'])
    trans = np.array(tab['trans'])

    return trans, ws



def getNormTrans(l, band='r', survey='sdss', bounds_error=False):
    """ 
    return the value of the normalized transmission function at wavelength l. 
    the function is from getNormTransFunc. 
    it is linear intrapolated to wavelength

    Params
    ------
    l (float): wavelength
    band='r'
    survey='sdss'
    bounds_error=False
            see interp1d
            If True, a ValueError is raised any time interpolation is attempted on
            a value outside of the range of x (where extrapolation is
            necessary). If False, out of bounds values are assigned `fill_value`.
            By default, an error is raised unless `fill_value="extrapolate"`.


    Return
    ------
    t (float): 
        value of the normalized transmission funciton
    """
    trans, ws = getNormTransFunc(band=band, survey=survey)

    f = interp1d(ws, trans, kind='linear', bounds_error=bounds_error, fill_value=0.)

    return max(f(l), 0)
