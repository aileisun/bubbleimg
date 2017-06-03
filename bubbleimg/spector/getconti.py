# getconti.py
# 12/05/2016 ALS 

"""
tools to isolate continuum 
"""

import numpy as np
import copy
from scipy.ndimage.filters import median_filter, generic_filter

from ..filters import getllambda


def decompose_conti_line_t2AGN(spec, lcoord, z):
    """
    decompose the spectrum of type 2 AGN into two components: continumm and emission line, by making the lines and running medium filter. This method is slower than a vanilla version medium filter but it estiamtes the continuum around the strong lines better. 

    See selectcont() for the line selection. 

    Params
    ------
    spec (array): spectrum, may have unit
    lcoord (array)
    z (float)
    toplot=False (bool)

    Return
    ------
    speccon
    specline
    lcoord
    """
    # main
    selcon = selectcont(spec, lcoord, z, AGN_TYPE=2, NLcutwidth=80., BLcutwidth=180., vacuum=True)

    speccon_nan = copy.copy(spec)
    speccon_nan[~selcon] = np.nan

    speccon = generic_filter(speccon_nan, np.nanmedian, size=300)

    specline = spec - speccon
    specline[selcon] = 0.


    if hasunit(spec):
        try:
            speccon.unit = spec.unit
        except:
            speccon = speccon * spec.unit
        try:
            specline.unit = spec.unit
        except:
            specline = specline * spec.unit

        units = [spec.unit, speccon.unit, specline.unit]

        if len(set(units)) != 1:
            raise Exception("[getconti] units are not identical")

    return speccon, specline, lcoord


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
    NL = np.array(['OII', 'OIII', 'OI', 'NII', 'NeV', 'SII', 'Ha', 'Hb', 'Hg', 'Hd', 'ArIII', 'NeIII', 'Dip', 'Mg', 'Na', 'He', 'HeI', 'HeII', 'FeVII'])
    llambdas = np.array([])
    for i in range(len(NL)):
        llambdas = np.append(llambdas, getllambda.getllambda(NL[i], vacuum=vacuum))
    llambdas = llambdas*(1.+z)
    selcut = (spec == 0)
    for i in range(llambdas.size):
        selline = np.absolute(xcoord-llambdas[i])<NLcutwidth/2.
        selcut = np.any([selcut, selline], axis=0)
    selkeep = np.logical_not(selcut)

    if AGN_TYPE == 1:
        BL = np.array(['Ha', 'Hb', 'Hg', 'Hd'])
        llambdas = np.array([])
        for i in range(len(BL)):
            llambdas = np.append(llambdas, getllambda.getllambda(BL[i], vacuum=vacuum))
        llambdas = llambdas*(1.+z)
        selcut = (spec==0)
        for i in range(llambdas.size):
            selline = np.absolute(xcoord-llambdas[i]) < BLcutwidth/2.
            selcut = np.any([selcut, selline], axis=0)
        selline = np.absolute(xcoord-getllambda.getllambda('Ha', vacuum=vacuum)*(1.+z))<200.
        selcut = np.any([selcut, selline], axis=0)
        selkeep = np.all([selkeep, np.logical_not(selcut)], axis=0)
    selkeep[np.all([np.arange(selkeep.size)>1870, np.arange(selkeep.size)<1905], axis=0)]=0
    return selkeep.astype('bool')


def getMedianFilteredConti_old(spec, lcoord, z):
    """
    PURPOSE: return median fitlered continuum, but in new lcoord system. 

    PARAMETERS: 
        spec (array)
        lcoord (array)
        z (float)
        toplot=False (bool)
    """

    # select continuum pixels
    selcon = selectcont(spec, lcoord, z, AGN_TYPE=2, NLcutwidth=80., BLcutwidth=180., vacuum=True)

    # median filter continuum
    speccon = median_filter(spec[selcon], size=30)
    lcoordcon = lcoord[selcon]

    # glue back the unit if spec has one
    try: 
        spec.unit
    except: 
        pass
    else: 
        speccon=speccon*spec.unit       

    return speccon, lcoordcon

