# inttools.py
# 06/02/2017 ALS 

"""
tool for integration
"""

import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import copy

import astropy.units as u
import astropy.constants as const

def calc_Fnu_in_band_from_fl(fl, ws, trans, ws_trans, isnormed=False):
    """
    calcualte: 

        Fnu =        (int fnu(l) T(l) dlnl) / (int T(l) dlnl)

           or

        Fnu = (int l^2/c * fl(l) T(l) dlnl) / (int T(l) dlnl)
    
    where fl is spectrum, T is filter transmission, l is wavelength

    Params
    ------
    fl (spec)     in units similar to "erg cm-2 s-1 AA-1"
    ws           in units similar to "AA-1"
    trans          dimensionless
    ws_trans    in units similar to "AA-1"
    isnormed=False (bool):
        if true, then assumes the transmission function is normalized such that int{ trans dlnl } = 1

    """

    # sanity check: that units are correct
    checks =[
            (fl/u.Unit("erg cm-2 s-1 AA-1")).to(u.dimensionless_unscaled).unit == u.dimensionless_unscaled, 
            (u.Quantity(trans)).unit == u.dimensionless_unscaled, 
            (ws/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            (ws_trans/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            ]
    if not all(checks):
        print checks
        raise Exception("[inttools] units of spectrum or filter response function is wrong")

    # convert to quantity
    fl = u.Quantity(fl)
    ws = u.Quantity(ws)

    # calculation
    fnu = (fl * (ws**2) / const.c).to(u.Unit("erg s-1 cm-2 Hz-1"))

    return calc_Fnu_in_band_from_fnu(fnu, ws, trans, ws_trans, isnormed=isnormed)



def calc_Fnu_in_band_from_fnu(fnu, ws, trans, ws_trans, isnormed=False):
    """
    calcualte: 

        Fnu =        (int fnu(l) T(l) dlnl) / (int T(l) dlnl)
    
    where fnu is spectrum, T is filter transmission, l is wavelength

    Params
    ------
    fnu (spec) in units similar to "erg s-1 cm-2 Hz-1"
    ws         in units similar to "AA-1"
    trans          dimensionless
    ws_trans   in units similar to "AA-1"
    transisnormed=False (bool):
        if true, then assumes the transmission function is normalized such that int{ trans dlnl } = 1
        it saves calculation time

    """

    # sanity check: that units are correct
    checks =[
            (fnu/u.Unit("erg s-1 cm-2 Hz-1")).unit == u.dimensionless_unscaled, 
            (u.Quantity(trans)).unit == u.dimensionless_unscaled, 
            (ws/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            (ws_trans/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            ]
    if not all(checks):
        print checks
        raise Exception("[inttools] units of spectrum or filter response function is wrong")

    # convert to quantity
    fnu = u.Quantity(fnu)
    ws = u.Quantity(ws)
    trans = u.Quantity(trans)
    ws_trans = u.Quantity(ws_trans)

    numerator = int_arr_times_arr_over_dlnx(arr1=fnu, xs1=ws, arr2=trans, xs2=ws_trans)

    if not isnormed:
        denominator = int_arr_over_dlnx(arr=trans, xs=ws_trans)
    else: 
        denominator = 1.

    Fnu = (numerator/denominator).to(u.Unit("erg s-1 cm-2 Hz-1"))

    return Fnu


def int_f_over_dlnx(f, x0=1., x1=2.):
    """
    Params
    ------
    f (float function to integrate)
    l0: start
    l1: end

    Return
    ------
    result
    """

    f_over_l = lambda l: f(l)/l

    result = integrate.quad(f_over_l, x0, x1)

    return result[0]


def int_arr_over_dlnx(arr, xs):
    """
    Params
    ------
    arr (array)
    xs (array)
    """

    lnxs = np.log(np.array(xs))

    return np.trapz(arr, x=lnxs)


def int_arr_times_f_over_dlnx(arr, f, xs):

    arr2 = f(xs)

    return int_arr_over_dlnx(arr*arr2, xs)


def int_arr_times_arr_over_dlnx(arr1, xs1, arr2, xs2, toexception=False):
    """
    integrate over arr1*arr2 dlnx, where xs1 and xs2 may be on different grids. 
    The integration is over d ln (xs1). arr2 is interpolated onto xs1 grids. 
    If xs1 is more extended than xs2 then arr2 outside of xs2 coverage is assumed to be 0. 

    If xs2 is more extended than xs1, raise error

    """
    arr1 = copy.copy(arr1)
    xs1 = copy.copy(xs1)

    if not arr1_contains_arr2(xs1, xs2) and toexception:
        raise Exception("[inttools] arr1 does not cover the domain of arr2")

    f = interp1d(xs2, arr2, kind='linear', bounds_error=False, fill_value=0.)
    return int_arr_times_f_over_dlnx(arr1, f, xs1)


def arr1_contains_arr2(arr1, arr2):

    return (max(arr1) >= max(arr2)) and (min(arr1) <= min(arr2))

