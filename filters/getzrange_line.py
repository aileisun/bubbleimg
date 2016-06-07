# getzrange_line.py
# ALS 2016/05/24

"""
find zrange for lines to be in/out the band. 
"""
import numpy as np
from astropy.table import Table, Row

def findzrange_wline_OIIIs(threshold=0.6):
    """
    Find the right red shift range such that both [OIII] 4960 and 5008 are within 
    throughput > 0.6*max range of r band. 
    """
    fileprefix='OIII'
    l1=getllambda(ion='OIII',lid=4960,vacuum=True) 
    l2=getllambda(ion='OIII',lid=5008,vacuum=True) 

    lmin=min(l1,l2)
    lmax=max(l1,l2)
    findzrange_line(fileprefix,lmin,lmax,threshold=threshold)



def findzrange_nline_HaNII(threshold=0.2):
    """
    Find the right red shift range such that both [OIII] 4960 and 5008 are within 
    throughput > 0.6*max range of r band. 
    """
    fileprefix='HaNII'
    l1=getllambda(ion='Ha')[0]
    l2,l3=getllambda(ion='NII')

    lmin=min(l1,l2,l3)
    lmax=max(l1,l2,l3)
    findzrange_line(fileprefix,lmin,lmax,threshold=threshold,inside=False)



def findzrange_nline_HaNIISII(threshold=0.2):
    """
    Find the right red shift range such that both [OIII] 4960 and 5008 are within 
    throughput > 0.6*max range of r band. 
    """
    fileprefix='HaNIISII'
    l1=getllambda(ion='Ha')[0]
    l2,l3=getllambda(ion='NII')
    l4,l5=getllambda(ion='SII')[1:]

    lmin=min(l1, l2, l3, l4, l5)
    lmax=max(l1, l2, l3, l4, l5)
    findzrange_line(fileprefix,lmin,lmax,threshold=threshold,inside=False)


def findzrange_line(linelist, l0, l1, inside=True, threshold=0.2):
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

    # setups
    if inside: 
        fileout='zrange_wline_'+linelist+'_'+'%.1f'%threshold+'.txt'
    else:
        fileout='zrange_nline_'+linelist+'_'+'%.1f'%threshold+'.txt'

    lmin, lmax = np.sort(np.array([l0,l1]))

    # make table
    tabout=Table([[],[],[],],names=('band','z0','z1'),dtype=('string','float','float'))

    for band in ['u','g','r','i','z']:

        # get filter boundary
        fb=Table.read('filterboundary_'+str(threshold)+'.txt',format='ascii')
        w1,w2=fb[fb['band']==band]['w1','w2'][0]

        # calculate corresponding redshift
        if inside:
            z0=w1/lmin-1.
            z1=w2/lmax-1.
        else: 
            z0=w2/lmin-1.
            z1=w1/lmax-1.

        tabout.add_row([band,z0,z1])

    # output
    tabout.write(fileout,format='ascii.fixed_width',delimiter='')
    return tabout


def getllambda(ion='OIII',lid=0,vacuum=True):
    """
    PURPOSE: Get the vaccum or air wavelength of a line in Angstrom. 

    INPUT:  getllambda(ion='OIII',lid=0,vacuum=True)
        ion: e.g., ion='OIII'
        lid:  e.g., lid=5008. If lid not specified then all lines of that ion will be returned. 
        vacuum: True if want vacuum wavelength, false if want air. 

    OUTPUT: line wavelengths in Angstrom
    """
    lid=int(lid)
    # e.g. ion='OIII', lid=5008
    linelist=Table.read('linelist.txt',format='ascii',delimiter='\t')
    if lid ==0:
        lam=linelist[linelist['Line']==ion]['Wavelength'].data
    else:
        if (linelist['Line']==ion).sum()==1:
            lam=linelist[linelist['Line']==ion]['Wavelength'].data
        elif (linelist['Line']==ion).sum()>1:
            lam=linelist[np.all([linelist['Line']==ion,linelist['Identifier']==lid],axis=0)]['Wavelength'].data
        else:
            raise NameError("[llambda] no line found")

    if vacuum:
        return lam
    else:
        from PyAstronomy import pyasl
        return pyasl.vactoair2(lam)


