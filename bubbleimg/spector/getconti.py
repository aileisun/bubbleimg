# getconti.py
# 12/05/2016 ALS 

"""
tools to isolate continuum 
"""

import numpy as np
import copy

from scipy.ndimage.filters import median_filter, generic_filter

from ..filters import getllambda
import linelist

def decompose_cont_line_t2AGN(spec, ls, z):
    """
    decompose the spectrum of type 2 AGN into two components: continumm and emission line, by making the lines and running medium filter. This method is slower than a vanilla version medium filter but it estiamtes the continuum around the strong lines better. 

    See selectcont() for the line selection. 

    Params
    ------
    spec (array): spectrum, may have unit
    ls (array)
    z (float)
    toplot=False (bool)

    Return
    ------
    speccon
    specline
    ls
    """
    # main
    selcon = selectcont(spec, ls, z, AGN_TYPE=2, NLcutwidth=80., BLcutwidth=180., vacuum=True)

    speccon_nan = copy.copy(spec)
    speccon_nan[~selcon] = np.nan

    speccon = generic_filter(speccon_nan, np.nanmedian, size=300)

    specline = spec - speccon
    specline[selcon] = 0.

    specline = inhereit_unit(specline, spec)
    speccon = inhereit_unit(speccon, spec)

    return speccon, specline, ls


def inhereit_unit(y, x):
    if hasunit(x):
        try:
            y.unit = x.unit
        except:
            y = y * x.unit
    return y


def hasunit(x):
    try:
        x.unit
    except:
        return False
    else:
        return True


def selectcont(spec, xcoord, z, AGN_TYPE=2, NLcutwidth=70., BLcutwidth=180., vacuum=True):
    """
    PURPOSE: select continuum w pixels using 1d spec
    PARAMETERS: 
        spec
        xcoord     [AA] wavelengths
        z
        AGN_TYPE=2
        NLcutwidth=70.
        BLcutwidth=180.
        vacuum=True
    """
    NL = linelist.narrowline
    BL = linelist.broadline

    llambdas = np.array([])
    for i in range(len(NL)):
        llambdas = np.append(llambdas, getllambda(NL[i], vacuum=vacuum))
    llambdas = llambdas*(1.+z)
    selcut = (spec == 0)
    for i in range(llambdas.size):
        selline = np.absolute(xcoord-llambdas[i])<NLcutwidth/2.
        selcut = np.any([selcut, selline], axis=0)
    selkeep = np.logical_not(selcut)

    if AGN_TYPE == 1:
        llambdas = np.array([])
        for i in range(len(BL)):
            llambdas = np.append(llambdas, getllambda(BL[i], vacuum=vacuum))
        llambdas = llambdas*(1.+z)
        selcut = (spec==0)
        for i in range(llambdas.size):
            selline = np.absolute(xcoord-llambdas[i]) < BLcutwidth/2.
            selcut = np.any([selcut, selline], axis=0)
        selline = np.absolute(xcoord-getllambda('Ha', vacuum=vacuum)*(1.+z)) < 200.
        selcut = np.any([selcut, selline], axis=0)
        selkeep = np.all([selkeep, np.logical_not(selcut)], axis=0)
    selkeep[np.all([np.arange(selkeep.size) > 1870, np.arange(selkeep.size) < 1905], axis=0)] = 0
    return selkeep.astype('bool')


