# getconti.py
# 12/05/2016 ALS 

"""
tools to isolate continuum 
"""

import numpy as np
import copy
import scipy.ndimage.filters as scif
import modelBC03

from ..filters import getllambda
from . import linelist

def decompose_cont_line_t2AGN(spec, ws, z, method='modelBC03'):
    """
    decompose the spectrum of type 2 AGN into two components: continumm and emission line. There are two methods: 

    method 'modelBC03':
        fit BC03 stellar population synthesis model to line masked continuum and use the best fit as the continuum model. 

    method 'running_median':
        mask the lines and do running medium filter. 

    See selectcont() for the line selection. 

    Params
    ------
    spec (array): spectrum, may have unit
    ws (array)
    z (float)
    toplot=False (bool)
    method = 'modelBC03'

    Return
    ------
    selcon
    speccon
    specline
    ws
    model
    """
    # main
    selcon = selectcont(spec, ws, z, AGN_TYPE=2, NLcutwidth=80., BLcutwidth=180., vacuum=True)

    if method == 'modelBC03':
        m = modelBC03.modelBC03(extinction_law='none')
        m.fit(ws=ws, spec=spec, z=z)
        speccon = m.bestfit_regrid
        model = m

    elif method == 'running_median':
        speccon = getcont_medianfilter(spec, selcon)
        model = None

    else:
        raise Exception("method is not recognized")

    specline = spec - speccon
    specline[selcon] = 0.

    specline = inherit_unit(specline, spec)
    speccon = inherit_unit(speccon, spec)

    return selcon, speccon, specline, ws, model


def getcont_medianfilter(spec, selcon):
    speccon_nan = copy.copy(spec)
    speccon_nan[~selcon] = np.nan
    speccon = scif.generic_filter(speccon_nan, np.nanmedian, size=300)

    return speccon


def inherit_unit(y, x):
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
