# getcontspec.py
# 12/05/2016 ALS 


"""
To get the line free continuum spectrum by masking the known lines and apply a medium band filter

see /Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/physical/contsubi.py
"""

import numpy as np
from scipy.ndimage.filters import median_filter

from ..filters import getllambda

# import sys
# import os
# sys.path.append(os.path.abspath('../filters'))
# import getllambda

def getMedianFilteredConti(spec, lcoord, z, toplot=False):
    """
    PURPOSE: return median fitlered continuum

    PARAMETERS: 
        spec (array)
        lcoord (array)
        z (float)
        toplot=False (bool)
    """

    # select continuum pixels
    selcon = selectcont(spec,lcoord,z,AGN_TYPE=2,NLcutwidth=80.,BLcutwidth=180.,vacuum=True)

    # median filter continuum
    specconmfiltered=median_filter(spec[selcon],size=30)
    lcoordcon=lcoord[selcon]

    if toplot:
        plt.clf()
        plt.plot(lcoord, spec, color='black')
        plt.plot(lcoordcon, specconmfiltered, color='r')

    # glue back the unit if spec has one
    try: spec.unit
    except: pass
    else: specconmfiltered=specconmfiltered*spec.unit       

    return specconmfiltered, lcoordcon



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
        selcut = np.any([selcut, selline],axis=0)
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
            selcut = np.any([selcut, selline],axis=0)
        selline = np.absolute(xcoord-getllambda.getllambda('Ha',vacuum=vacuum)*(1.+z))<200.
        selcut = np.any([selcut, selline],axis=0)
        selkeep = np.all([selkeep, np.logical_not(selcut)], axis=0)
    selkeep[np.all([np.arange(selkeep.size)>1870, np.arange(selkeep.size)<1905], axis=0)]=0
    return selkeep.astype('bool')

