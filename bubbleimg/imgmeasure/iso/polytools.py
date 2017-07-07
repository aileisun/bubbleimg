# py
# ALS 2016/04/29

"""
Contain functions for shape measurements on polygons. 

Does not support astropy unit

References: 
https://en.wikipedia.org/wiki/Shape_factor_(image_analysis_and_microscopy)
https://en.wikipedia.org/wiki/Feret_diameter
PEER REVIEW Particle Shape Factors and Their Use in Image Analysis Part 1: Theory by Eric Olson


Notation
--------
    polys: list of poly
    poly: (N*2) ndarray of double
        list of corner coordinates (r,c)=(y,x) of the polygon.


Next developement steps: 
    1. build FeretD_max.
"""

import numpy as np
import skimage.measure as skmeasure
import collections

import astropy.table as at

def ShapeParamsTab_from_contours(contours, xc, yc):
    """
    return a table of measurements from input contours (everything in pix unit) e.g. :

    area_pix    dmax_pix      rmax_pix      dper_pix     theta_dmax    theta_rmax    theta_dper   aspectr
    -------- ------------- ------------- ------------- ------------- ------------- -------------- -------
         8.0 3.60555127546 1.80277563773 3.60555127546 33.6900668147 146.309932696 -56.3099331853     1.0

    """
    cols = ['area_pix', 'dmax_pix', 'rmax_pix', 'dper_pix', 'theta_dmax', 'theta_rmax', 'theta_dper', 'aspectr']
    tab = at.Table(names=cols)

    if len(contours)>0:
        area = NetPolygonsArea(contours)
        dferetmax, theta_dferetmax = FeretD_max(contours)
        rferetmax, theta_rferetmax = FeretR_max(contours, xc, yc)
        dferetper, theta_dferetper = FeretD_max90(contours)
        aspectr = FeretAspectRatio(contours)
        row = [area, dferetmax, rferetmax, dferetper, theta_dferetmax, theta_rferetmax, theta_dferetper, aspectr]
    else: 
        row = [0., 0., 0., 0., np.nan, np.nan, np.nan, np.nan]

    tab.add_row(row)

    return tab

# def ShapeParamsdict_from_contours(contours, xc, yc):
#     """
#     return a dictionary of key:value pairs of ShapeParams measured from contours
#     """
#     cols = ['area', 'dferetmax', 'theta_dferetmax', 'rferetmax', 'theta_rferetmax', 'dferetper', 'theta_dferetper', 'aspectr']
#     if len(contours)>0:
#         area = NetPolygonsArea(contours)
#         dferetmax, theta_dferetmax = FeretD_max(contours)
#         rferetmax, theta_rferetmax = FeretR_max(contours, xc, yc)
#         dferetper, theta_dferetper = FeretD_max90(contours)
#         aspectr = FeretAspectRatio(contours)
#         values=[area, dferetmax, theta_dferetmax, rferetmax, theta_rferetmax, dferetper, theta_dferetper, aspectr]
#     else: 
#         values=[0., 0., np.nan, 0., np.nan, 0., np.nan, np.nan]

#     paramsdict = collections.OrderedDict(zip(cols, values))
#     return paramsdict


def find_centercontours(img, threshold, xc, yc, radius=3):
    contours = find_contours(img, threshold)   
    return select_center_contours(contours, xc, yc, radius=radius) 


def find_largecontours(img, threshold, minarea):
    """ Find all isophote patch larger than certain area. See select_large_contours for details.  """

    contours = find_contours(img, threshold)

    return select_large_contours(contours, minarea)


def find_contours(img, threshold, tocomplete=True):
    """
    Return all isophote contours. Each contours is (N*2) ndarray of double. 
    High contours are counter-clockwise points and low contours are clockwise 
    points.

    if tocomplete is True, even contours at edges are always closed. 
    """
    # negative edge padding to close all contours
    if tocomplete: # to complete contours
        wpad=2
        img_pad = np.pad(img, [[wpad, wpad], [wpad, wpad]], mode='constant', constant_values=-1.e300)
        contours = skmeasure.find_contours(img_pad, threshold, fully_connected='low', positive_orientation='high')
        contours = [contour-wpad for contour in contours]
    else: 
        contours = skmeasure.find_contours(img, threshold, fully_connected='low', positive_orientation='high')

    return contours


def select_large_contours(contours, minarea):
    """
    Find all patch that have areas above minarea. A patch is a highcontour plus all of its enclosing lowcontours. Here the holes are considered. 
    """
    large_highcontours = select_large_highcontours(contours, minarea)
    
    results = []

    for highcontour in large_highcontours: 
        patch = get_patches_of_highcontours(contours, [highcontour])
        if NetPolygonsArea(patch) >= minarea:
            results = results+patch

    return unique_polys(results)


def select_large_highcontours(contours, minarea):
    """
    Find all high contours that have areas above minarea (in pixel). The holes within the high contours are ignored. 
    """
    results = []

    for contour in contours:
        if SignedPolygonArea(contour) >= minarea:
            results = results+[contour]
    return unique_polys(results)


def select_center_contours(contours, xc, yc, radius=3):

    highcontours = select_center_highcontours(contours, xc, yc, radius=radius)
    allcontours = get_patches_of_highcontours(contours, highcontours)

    return allcontours


def select_center_highcontours(contours, xc, yc, radius=3):
    """
    select all the high contours in contours that overlaps with the center region, which is defeined as a square region (xc +/- radius, yc +/- radius). 

    Find all high contours that enclose the centroid xc, yc. +/- radius
    # if skmeasure.points_in_poly([[yc, xc]], contour):
    """
    # setting
    carea = np.pi * radius**2
    ptc = np.array([xc, yc])

    highcontours = select_highcontours(contours)

    ccontours = []

    for contour in highcontours:
        if SignedPolygonArea(contour) >= carea:
            if skmeasure.points_in_poly(np.array([ptc]), contour)[0]:
                ccontours = ccontours+[contour]
        else: 
            if contour_is_close_to_point(contour, ptc, distance=radius):
                ccontours = ccontours+[contour]

    return ccontours


def contour_is_close_to_point(contour, pt, distance=1.):

    dissq = ((pt[0]-contour[:, 0])**2 + (pt[1]-contour[:, 1])**2)

    result = np.any(dissq <= distance**2)
    return result



# def select_center_highcontours(contours, xc, yc, radius=3):
#     """
#     select all the high contours in contours that overlaps with the center region, which is defeined as a square region (xc +/- radius, yc +/- radius). 

#     Find all high contours that enclose the centroid xc, yc. +/- radius
#     # if skmeasure.points_in_poly([[yc, xc]], contour):
#     """
#     # define neighboring points
#     xns = xc + np.arange(-radius, radius+1, step=1)
#     yns = yc + np.arange(-radius, radius+1, step=1)
#     pns = [[yn, xn] for yn in yns for xn in xns] 

#     highcontours = select_highcontours(contours)
#     ccontours = []

#     for contour in highcontours:
#         if np.any(skmeasure.points_in_poly(pns, contour)):
#             ccontours = ccontours+[contour]

#     return ccontours


def select_highcontours(contours):
    """
    select high contours within contours
    """
    highcontours = []
    for contour in contours:
        if SignedPolygonArea(contour) > 0.:
            highcontours = highcontours+[contour]
    return highcontours


def select_lowcontours(contours):
    """
    select low contours within contours
    """
    lowcontours = []
    for contour in contours:
        if SignedPolygonArea(contour) < 0.:
            lowcontours = lowcontours+[contour]
    return lowcontours


def select_enclosedcontours(contours, outcontours):
    """
    Among 'contours', select the ones that are enclosed by 'outcontours'. The contours that are identical to any of the outcontours are not included. 
    """
    contours = to_polys(contours)
    outcontours = to_polys(outcontours)

    incontours = []
    for contour in contours:
        isinsides = np.array([isinside_polygon(contour, outcontour) for outcontour in outcontours])
        isequals = np.array([isequal_polygon(contour, outcontour) for outcontour in outcontours])
        if any(isinsides) and all(~isequals): 
            incontours = incontours+[contour]
    return incontours


def to_polys(poly):
    if ispolys(poly): 
        return poly
    else:
        if ispolys([poly]):
            return [poly]
        else:
            raise ValueError("input not polys or poly")


def get_patches_of_highcontours(contours, highcontours):
    """ find all the low contours enclosed in the highcontours, return highcontours + enclosed low contours """

    if not ispolys(highcontours):
        raise ValueError("[polytools] highcontours is not polys")
    if not ispolys(contours):
        raise ValueError("[polytools] contours is not polys")

    lowcontours = select_lowcontours(contours)
    lowcontours = select_enclosedcontours(lowcontours, highcontours)

    allcontours = highcontours+lowcontours
    allcontours = unique_polys(allcontours)

    return allcontours



def unique_polys(polys):
    """
    Return unique polys
    """
    uniqpolys = []
    for poly in polys:
        if not any([isequal_polygon(p, poly) for p in uniqpolys]):
            uniqpolys.append(poly)
    return uniqpolys


def ispolys(polys):
    """
    tell if 'polys' is a list of polygons, which are np arrays of (x,y) points.
    """
    if isinstance(polys, (list, np.ndarray)):
        if np.all([ispoly(poly) for poly in polys]):
            return True
        else:
            return False
    else:
        return False


def ispoly(poly):
    """
    return True if poly is a np array of shape (N,2)
    """
    if isinstance(poly, np.ndarray):
        if len(poly.shape) == 2 and poly.shape[-1] == 2:
            return True
        else:
            return False
    else:
        return False


def isequal_polygon(poly1, poly2):
    """ tell if these two polygons are equal """
    return np.array_equal(poly1, poly2)


def isinside_polygon(poly1, poly2):
    """ tell if poly1 is enclosed by poly2, i.e. all the points are inside.
    Identical polys are not inside each other """
    return np.all(skmeasure.points_in_poly(poly1, poly2))


def NetPolygonsArea(polys):
    """ 
    Calculate the net area of a groupd of polygons. Counter-clock wise 
    polygons have positive areas while clock-wise ones have negative. 

    Parameters
    ----------
    polys : list of polygons, each is a (N*2) array of corner coordinates
        list of polygons to use 

    Returns
    -------
    area : float
    """
    area = 0.0
    for i in range(len(polys)):
        area += SignedPolygonArea(polys[i])
    return area


def SignedPolygonArea(poly):
    """ 
    Calculate the signed area inside a polygon using Shoelace formula. 
    Positive for clock wise (on x-y plane) polygons. 

    Applied correction to take into account the fact that the four corners of the box can be dropped because of the interpolation. After the correction, we get exact area for the contours of a binary file where the threshold is set to 0.5.  See test_polytools.py

    Parameters
    ----------
    poly : (N*2) ndarray of double
        list of corner coordinates (r,c) = (y,x) of the polygon 
        N is the number of corners. 

    Returns
    -------
    area : float

    Reference
    -------
    https://en.wikipedia.org/wiki/Shoelace_formula
    http://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    """
    n = len(poly) # of poly
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += poly[i][0] * poly[j][1]
        area -= poly[j][0] * poly[i][1]
    area = area / 2.0

    if area > 0:
        cornercorrection = 0.5
    elif area < 0:
        cornercorrection = -0.5

    return area + cornercorrection


def FeretDiameter(sth,theta):
    """
    Return the Feret's diamter along the direction of theta (y of x.) of the object 'sth'. 'sth' can either be polys (list of poly) or poly (array of shape (N,2)). Only high contours are considered. 

    Parameters
    ----------
    sth: either polys or poly
        polys: list of poly
        poly: (N*2) ndarray of double
            list of corner coordinates (r,c) of the polygon.

    theta: (float)
        the angle of the direction being measured. y of x in degrees. 

    Returns
    ----------
    d_feret: (float)
        the Feret's diameter in pix
    """
    # define unit vector (r,c)=(y,x) along the direction of theta

    listpoints = np.ndarray([0,2])

    if ispolys(sth):
        polys = sth
        for poly in polys:
            if SignedPolygonArea(poly)>0.:
                listpoints = np.append(listpoints,poly,axis=0)
    elif ispoly(sth):
        poly=sth
        if SignedPolygonArea(poly)>0.:
            listpoints = np.append(listpoints,poly,axis=0)
    else:
        raise ValueError('input is neither a poly nor polys')

    unitvec=np.array([np.sin(np.radians(theta)),np.cos(np.radians(theta))])
    d = listpoints.dot(unitvec)
    dotmax=np.max(d,axis=0)
    dotmin=np.min(d,axis=0)
    return np.absolute(dotmax-dotmin)


def FeretRadius(sth, theta, xc, yc):
    """
    Like FeretDiameter but return the max radius (distance between edge and centroid) along the direction theta. 
    Parameters
    ----------
    sth: either polys or poly

    theta: (float)
        the angle of the direction being measured. y of x in degrees. 

    xc, yc: (float)
        the centroid coordinate

    Returns
    ----------
    r_feret: (float)
        the Feret's diameter in pix
    """
    # define unit vector (r,c)=(y,x) along the direction of theta

    listpoints=np.ndarray([0,2])

    if ispolys(sth):
        polys=sth
        for poly in polys:
            if SignedPolygonArea(poly)>0.:
                listpoints=np.append(listpoints,poly,axis=0)
    elif ispoly(sth):
        poly=sth
        if SignedPolygonArea(poly)>0.:
            listpoints=np.append(listpoints,poly,axis=0)
    else:
        raise ValueError('input is neither a poly nor polys')

    unitvec=np.array([np.sin(np.radians(theta)),np.cos(np.radians(theta))])
    pc=np.array([yc,xc])
    rs=(listpoints-pc).dot(unitvec)
    return np.max(rs,axis=0)


def FeretR_max(sth, xc, yc):
    """
    Return the maximum Feret's radius and its orientaion (y of x). 

    Parameters
    ----------
    sth: either polys or poly
        polys: list of poly
        poly: (N*2) ndarray of double
            list of corner coordinates (r,c) of the polygon.

    xc, yc: image center

    Returns
    ----------
    d_feret: (float)
        the Feret's diameter in pix
    theta: (float)
        the angle of the direction being measured. y of x in degrees. 
    """
    from scipy.optimize import minimize

    # init
    thetas=np.linspace(0., 360., num=15, endpoint=False)
    rfs=FeretRadius(sth ,thetas, xc, yc)
    theta_init=thetas[np.argmax(rfs)]
    tomin = lambda theta, sth: -FeretRadius(sth , theta, xc, yc)

    # fit
    results=minimize(tomin, theta_init,method='Powell', args=(sth,))
    if not results['success']: raise NameError('minimization failed')

    # return
    theta_max=results['x'].sum()
    rf_max=FeretRadius(sth, theta_max, xc, yc)
    return rf_max, theta_max


def FeretD_max(sth):
    """
    Return the maximum Feret's diameter and its orientaion (y of x). 

    Parameters
    ----------
    sth: either polys or poly
        polys: list of poly
        poly: (N*2) ndarray of double
            list of corner coordinates (r,c) of the polygon.

    Returns
    ----------
    d_feret: (float)
        the Feret's diameter in pix
    theta: (float)
        the angle of the direction being measured. y of x in degrees. 
    """
    from scipy.optimize import minimize

    # init
    thetas = np.linspace(0., 180., num = 15, endpoint=False)
    dfs = FeretDiameter(sth,thetas)
    theta_init=thetas[np.argmax(dfs)]
    tomin = lambda theta, sth: -FeretDiameter(sth, theta)

    # fit
    results = minimize(tomin, theta_init,method='Powell', args=(sth,))
    if not results['success']: raise NameError('minimization failed')

    # return
    theta_max = results['x'].sum()
    df_max = FeretDiameter(sth,theta_max)
    return df_max, theta_max


def FeretD_min(sth):
    """
    Return the minimum Feret's diameter and its orientaion (y of x). 

    Parameters
    ----------
    sth: either polys or poly
        polys: list of poly
        poly: (N*2) ndarray of double
            list of corner coordinates (r,c) of the polygon.

    Returns
    ----------
    d_feret: (float)
        the Feret's diameter in pix
    theta: (float)
        the angle of the direction being measured. y of x in degrees. 
    """
    from scipy.optimize import minimize

    # init
    thetas = np.linspace(0., 180., num=120, endpoint=False)
    dfs = FeretDiameter(sth, thetas)
    theta_init = thetas[np.argmin(dfs)]
    tomin = lambda theta, sth: FeretDiameter(sth, theta)

    # fit
    results = minimize(tomin, theta_init, method='Powell', args=(sth,))

    if not results['success']:
        raise NameError('minimization failed')

    # return
    theta_min = results['x'].sum()
    df_min = FeretDiameter(sth, theta_min)
    return df_min, theta_min


def FeretD_max90(sth):
    """
    Return Feret's diameter along the direction perpendicular to the maximum
    Feret's diameter.
    """
    df_max, theta_max = FeretD_max(sth)
    theta_max90 = theta_max-90.
    df_max90 = FeretDiameter(sth, theta_max90)
    return df_max90, theta_max90


def FeretAspectRatio(sth):
    """
    Return the aspect ratio as defined by FeretD_max90/FeretD_max
    where FeretD_max is the maximum Feret's Diamter and FeretD_max90 is the
    Feret's diameter perpendicular to the maximum's direction.
    """
    return FeretD_max90(sth)[0]/FeretD_max(sth)[0]


def mom_ellipseparams_poly(poly):
    """
    Return the (estimated) moment ellipse params of a polygon.

    Parameters
    ----------
    poly : (N*2) ndarray of double
        list of corner coordinates (r,c) of the polygon.

    Returns
    ----------
    params: list
        [m0, xc, yc, a, b, theta]
        where m0 is the total area
    """
    from copy import deepcopy
    from skimage.measure import points_in_poly
    import moment
    poly = deepcopy(poly)

    xmin, xmax = np.min(poly[:, 1], axis=0.), np.max(poly[:, 1], axis=0.)
    ymin, ymax = np.min(poly[:, 0], axis=0.), np.max(poly[:, 0], axis=0.)

    poly[:, 1] = poly[:, 1]-xmin
    poly[:, 0] = poly[:, 0]-ymin

    img = np.zeros([int(ymax-ymin), int(xmax-xmin)])

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if points_in_poly([[i,j]],poly):
                img[i,j]=1.

    if np.sum(img)>0.:
        params=moment.mom_ellipseparams(img)
        params[1]=xmin+params[1]
        params[2]=ymin+params[2]
        return params
    else: 
        return [0., 0., 0., 0., 0., 0.]



# def select_center_highcontours(contours, xc, yc, radius=3, areallimit_ctr=0.):
#     """
#     Find all high contours that enclose the centroid xc, yc. +/- radius
#     # if skmeasure.points_in_poly([[yc, xc]], contour):
#     """
#     # define neighboring points
#     xns = xc + np.arange(-radius, radius+1, step=1)
#     yns = yc + np.arange(-radius, radius+1, step=1)
#     pns = [[yn, xn] for yn in yns for xn in xns] 

#     ccontours = []
#     for contour in contours:
#         if SignedPolygonArea(contour) > radius**2:
#             if np.any(skmeasure.points_in_poly(pns, contour)):
#                 ccontours = ccontours+[contour]
#         elif SignedPolygonArea(contour) > areallimit_ctr:
#             if np.absolute(contour[:, 0]-yc).min() < radius: 
#                 if np.absolute(contour[:, 0]-xc).min() < radius: 
#                     if np.any(skmeasure.points_in_poly(pns, contour)):
#                         ccontours = ccontours+[contour]
#     return ccontours

# def find_realisocontours(img, threshold, minarea, xc, yc):
#     """
#     Return isophote high contours that are either centered on the image
#     center or has area larger than minarea, as well as their enclosed low
#     contours.
#     """
#     contours = find_contours(img, threshold)

#     ccontours = select_center_highcontours(contours, xc, yc)
#     largecontours = select_large_highcontours(contours, minarea)
#     highrealcontours = ccontours+largecont    ours

#     allcontours = get_patch_of_highcontour(contours, highcontours)

#     return allcontours
