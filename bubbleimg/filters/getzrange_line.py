# getzrange_line.py
# ALS 2016/05/24

"""
find zrange for lines to be in/out the band. 
"""

import numpy as np
from astropy.table import Table, Row

from . import filtertools
from . import getllambda


def findzrange_wline_OIIIs(threshold=0.6, survey='sdss'):
    """
    Find the right red shift range such that both [OIII] 4960 and 5008 are within 
    throughput > 0.6*max range of r band. 
    """
    fileprefix = 'OIII'
    # l1 = getllambda.getllambda(ion='OIII',lid=4960,vacuum=True) 
    l2 = getllambda.getllambda(ion='OIII',lid=5008,vacuum=True) 

    lmin = l2
    lmax = l2
    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, survey=survey)



def findzrange_wline_HaNII(threshold=0.6, survey='sdss'):
    """
    Find the right red shift range such that all of Ha + 2 NII are in the band
    """
    fileprefix = 'HaNII'
    l1 = getllambda.getllambda(ion='Ha')[0]
    l2 = getllambda.getllambda(ion='NII', lid=6550)
    l3 = getllambda.getllambda(ion='NII', lid=6585)

    lmin = min(l1, l2, l3)
    lmax = max(l1, l2, l3)
    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, survey=survey)



def findzrange_nline_HaNII(threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that none of Ha + 2 NII are in the band
    """
    fileprefix = 'HaNII'
    l1 = getllambda.getllambda(ion='Ha')[0]
    l2 = getllambda.getllambda(ion='NII', lid=6550)
    l3 = getllambda.getllambda(ion='NII', lid=6585)

    lmin = min(l1,l2,l3)
    lmax = max(l1,l2,l3)
    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, inside=False, survey=survey)


def findzrange_nline_HaNIISII(threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that Halpha NII an SII are outside specified bands
    """
    fileprefix = 'HaNIISII'
    l1 = getllambda.getllambda(ion='Ha')[0]
    l2 = getllambda.getllambda(ion='NII', lid=6550)
    l3 = getllambda.getllambda(ion='NII', lid=6585)
    l4,l5 = getllambda.getllambda(ion='SII')[1:]

    lmin = min(l1, l2, l3, l4, l5)
    lmax = max(l1, l2, l3, l4, l5)
    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, inside=False,  survey=survey)


def findzrange_nline_OII(threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that the OII 3727 are outside specified bands
    """
    fileprefix = 'OII'

    l1 = getllambda.getllambda(ion='OII')[0]
    findzrange_line(fileprefix, l1, l1, threshold=threshold, inside=False,  survey=survey)


def findzrange_nline_OIINeIII(threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that the NeIII 3869.81 are outside specified bands
    """
    fileprefix = 'OIINeIII'

    l1 = getllambda.getllambda(ion='OII')[0]
    l2 = getllambda.getllambda(ion='NeIII')[0]

    lmin = min(l1, l2)
    lmax = max(l1, l2)

    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, inside=False,  survey=survey)


def findzrange_nline_OIII(threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that OIII 4959 and 5007 are outside specified bands
    """

    fileprefix = 'OIII'
    l1 = getllambda.getllambda(ion='OIII',lid=4960,vacuum=True) 
    l2 = getllambda.getllambda(ion='OIII',lid=5008,vacuum=True) 

    lmin = min(l1, l2)
    lmax = max(l1, l2)

    findzrange_line(fileprefix, lmin, lmax, threshold=threshold, inside=False,  survey=survey)


def findzrange_line(linelist, l0, l1, inside=True, threshold=0.2, survey='sdss'):
    """
    Find the right red shift range such that lines are inside/outside of
    throughput > treshold*max range of band. 

    Note that if choose outside (inside=False), in the file z0>z1 for later 
    convenience. 


    Parameters
    -----
    linelist: string
        like OIII or HaNIISII

    l0, l1: float 
        line wavelength range in angstrom, doesn't need to be in order

    inside: bool
        if True, then return the z range that contains line in the band, 
        otherwise the range that have non of the lines in the band. 
    """

    # setup file paths
    localpath = filtertools.getlocalpath()
    filefb = localpath+survey+'/'+'filterboundary_'+str(threshold)+'.txt'
    if inside: 
        fileout = localpath+survey+'/'+'zrange_wline_'+linelist+'_'+'%.1f'%threshold+'.txt'
    else:
        fileout = localpath+survey+'/'+'zrange_nline_'+linelist+'_'+'%.1f'%threshold+'.txt'

    # setup params
    bands = filtertools.surveybands[survey]

    # make table
    lmin, lmax = np.sort(np.array([l0,l1]))

    tabout = Table([[],[],[],], names=('band','z0','z1'), dtype=('S1', 'f8', 'f8'))

    for band in bands:

        # get filter boundary
        fb = Table.read(filefb,format='ascii')
        w1,w2 = fb[fb['band']==band]['w1','w2'][0]

        # calculate corresponding redshift
        if inside:
            z0 = w1/lmin-1.
            z1 = w2/lmax-1.
        else: 
            z0 = w2/lmin-1.
            z1 = w1/lmax-1.

        tabout.add_row([band,z0,z1])

    # output
    tabout.write(fileout, format='ascii.fixed_width', delimiter='', overwrite=True)
    return tabout


    