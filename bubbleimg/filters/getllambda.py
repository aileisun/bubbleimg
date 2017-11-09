# getllambda.py
# ALS 2016/12/05

"""
Get wavelength of line in AA (in vacuum by default). The data is read from file linelist.txt. 
"""

import numpy as np
import re

import astropy.table as at
from PyAstronomy import pyasl

from . import filtertools


def getllambda(ion='OIII', lid=0, vacuum=True):
    """
    PURPOSE: Get the vaccum or air wavelength of a line in Angstrom. 

    getllambda(ion='OIII',lid=0,vacuum=True)
    
    Params
    ------
    ion (str): 
        e.g., ion='OIII'
            \or 
        e.g., ion='OIII5008' and not specify lid
            if lid is not specified then use the trialing numbers as lid

    lid (int): 
        e.g., lid=5008. :
        If lid not specified then all lines of that ion will be returned. 

    vacuum: 
        True if want vacuum wavelength, false if want air. 

    Return
    ------
    line wavelengths in Angstrom
    """
    localpath = filtertools.getlocalpath()
    filein = localpath+'linelist.txt'
    linelist = at.Table.read(filein,format='ascii',delimiter='\t')

    lid=int(lid)

    # if ion has trailing wavelength number then parse it to lid
    m = re.search(r'\d+$', ion)
    if m is not None: # there is trailing number
        if lid == 0:
            lid = int(m.group())
            ion = str(ion[:-len(str(lid))])
        else: 
            raise Exception("lid specified twice")

    if lid == 0:
        lam=linelist[linelist['Line']==ion]['Wavelength'].data
    else:
        if (linelist['Line']==ion).sum()==1:
            lam=linelist[linelist['Line']==ion]['Wavelength'].data
        elif (linelist['Line']==ion).sum()>1:
            lam=linelist[np.all([linelist['Line']==ion,linelist['Identifier']==lid],axis=0)]['Wavelength'].data
        else:
            raise NameError("[llambda] no line found for {} {}".format(ion, str(lid)))

    if vacuum:
        return lam
    else:
        return pyasl.vactoair2(lam)
