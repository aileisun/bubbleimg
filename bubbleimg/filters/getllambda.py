# getllambda.py
# ALS 2016/12/05

"""
Get wavelength of line in AA (in vacuum by default). The data is read from file linelist.txt. 
"""

import numpy as np

import astropy.table as at
from PyAstronomy import pyasl

import filtertools


def getllambda(ion='OIII',lid=0,vacuum=True):
    """
    PURPOSE: Get the vaccum or air wavelength of a line in Angstrom. 

    INPUT:  getllambda(ion='OIII',lid=0,vacuum=True)
        ion: e.g., ion='OIII'
        lid:  e.g., lid=5008. If lid not specified then all lines of that ion will be returned. 
        vacuum: True if want vacuum wavelength, false if want air. 

    OUTPUT: line wavelengths in Angstrom
    """
    localpath = filtertools.getlocalpath()
    filein = localpath+'linelist.txt'

    lid=int(lid)
    # e.g. ion='OIII', lid=5008
    linelist = at.Table.read(filein,format='ascii',delimiter='\t')
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
        return pyasl.vactoair2(lam)


