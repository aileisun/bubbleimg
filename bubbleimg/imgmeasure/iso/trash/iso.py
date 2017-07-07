# iso.py
# ALS 2015/04/22

"""
Make and return isophotal measurements
"""

import numpy as np

import polytools
from .. import standards



def iso_ShapeParams(img, threshold, areallimit):
    """
    Make shape measurements to the iso contours. The contour at the center 
    and the contours larger than an area of areallimit are considered. 

    Parameters
    ------
    img (2d np array)
    threshold (float): threshold of the isophote

    Returns
    ------
    paramsdict: dictionary
    ---
    area:            net area in the centroid "isophote" (minus the holes). 
    dferetmax:       maximum feret's diameter
    theta_dferetmax: orientation of the dferetmax in degree y of x
    rferetmax:       maximum radius measured from center
    theta_rferetmax: orientation of the rferetmax in degree y of x
    dferetper:       feret's diameter perpendicular to dferetmax
    theta_dferetper: orientation of the dferetper in degree y of x
    aspectr:         aspect ratio
    
    DESCRIPTION
    ------
    see polytools. ShapeParamsdict_from_contours.
    """
    xc, yc = standards.get_img_xycenter(img)

    # find contour
    contours = polytools.find_realisocontours(img, threshold, areallimit, xc, yc)
    # calculate
    paramsdict = polytools.ShapeParamsdict_from_contours(contours, xc, yc)
    return paramsdict



