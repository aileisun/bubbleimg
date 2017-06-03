# inttools.py
# 06/02/2017 ALS 

"""
tool for integration
"""

import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.constants as const

def calc_Fnu_in_band_from_fl(fl, ls, trans, ls_trans):
    """
    calcualte: 

        Fnu =        (int fnu(l) T(l) dlnl) / (int T(l) dlnl)

           or

        Fnu = (int l^2/c * fl(l) T(l) dlnl) / (int T(l) dlnl)
    
    where fl is spectrum, T is filter transmission, l is wavelength

    Params
    ------
    fl (spec)     in units similar to "erg cm-2 s-1 AA-1"
    ls         in units similar to "AA-1"
    trans          dimensionless
    ls_trans   in units similar to "AA-1"

    """

    # sanity check: that units are correct
    checks =[
            (fl/u.Unit("erg cm-2 s-1 AA-1")).unit == u.dimensionless_unscaled, 
            (u.Quantity(trans)).unit == u.dimensionless_unscaled, 
            (ls/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            (ls_trans/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            ]
    if not all(checks):
        print checks
        raise Exception("[inttools] units of spectrum or filter response function is wrong")

    # convert to quantity
    fl = u.Quantity(fl)
    ls = u.Quantity(ls)

    # calculation
    fnu = (fl * (ls**2) / const.c).to(u.Unit("erg s-1 cm-2 Hz-1"))

    return calc_Fnu_in_band_from_fnu(fnu, ls, trans, ls_trans)



def calc_Fnu_in_band_from_fnu(fnu, ls, trans, ls_trans):
    """
    calcualte: 

        Fnu =        (int fnu(l) T(l) dlnl) / (int T(l) dlnl)
    
    where fnu is spectrum, T is filter transmission, l is wavelength

    Params
    ------
    fnu (spec) in units similar to "erg s-1 cm-2 Hz-1"
    ls         in units similar to "AA-1"
    trans          dimensionless
    ls_trans   in units similar to "AA-1"

    """

    # sanity check: that units are correct
    checks =[
            (fnu/u.Unit("erg s-1 cm-2 Hz-1")).unit == u.dimensionless_unscaled, 
            (u.Quantity(trans)).unit == u.dimensionless_unscaled, 
            (ls/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            (ls_trans/u.Unit("AA")).unit == u.dimensionless_unscaled, 
            ]
    if not all(checks):
        print checks
        raise Exception("[inttools] units of spectrum or filter response function is wrong")

    # convert to quantity
    fnu = u.Quantity(fnu)
    ls = u.Quantity(ls)
    trans = u.Quantity(trans)
    ls_trans = u.Quantity(ls_trans)

    numerator = int_arr_times_arr_over_dlnx(arr1=fnu, xs1=ls, arr2=trans, xs2=ls_trans)
    denominator = int_arr_over_dlnx(arr=trans, xs=ls_trans)

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


def int_arr_times_arr_over_dlnx(arr1, xs1, arr2, xs2):
    """
    integrate over arr1*arr2 dlnx, where their x grids may be different
    use arr1, xs1 as the integral reference, interpolate arr2 onto xs1
    """

    f = interp1d(xs2, arr2, kind='linear', bounds_error=False, fill_value=0.)

    return int_arr_times_f_over_dlnx(arr1, f, xs1)
