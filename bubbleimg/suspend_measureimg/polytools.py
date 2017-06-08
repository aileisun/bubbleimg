# polytools.py
# ALS 2016/04/29

"""
Contain functions for shape measurements on polygons. 

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

from .. import standards


def find_realisocontours(img, threshold, areallimit):
    """
    Return isophote high contours that are either centered on the image
    center or has area larger than areallimit, as well as their enclosed low
    contours.
    """
    # xc, yc = 0.5*(np.array(img.shape))
    xc, yc = standards.get_img_xycenter(img)

    contours = find_isocontours(img, threshold)

    ccontours = select_centercontours(contours, xc, yc)
    largecontours = select_largecontours(contours, areallimit)
    highrealcontours = ccontours+largecontours

    lowcontours = select_lowcontours(contours)
    lowrealcontours = select_enclosedcontours(lowcontours, highrealcontours)

    realcontours = highrealcontours+lowrealcontours
    realcontours = unique_polys(realcontours)
    return realcontours


def find_isocontours(img, threshold, tocomplete=True, nx=64, ny=64):
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


def unique_polys(polys):
    """
    Return unique polys
    """
    uniqpolys = []
    for poly in polys:
        if not any([isequal_polygon(p, poly) for p in uniqpolys]):
            uniqpolys.append(poly)
    return uniqpolys


def select_largecontours(contours, areallimit):
    """
    Find all high contours that have areas above areallimit (in pixel).
    """
    largecontours = []
    for contour in contours:
        if SignedPolygonArea(contour) > areallimit:
            largecontours = largecontours+[contour]
    return largecontours


def select_centercontours(contours, xc, yc, radius=3, areallimit_ctr=3):
    """
    Find all high contours that enclose the centroid xc, yc. +/- radius
    # if skmeasure.points_in_poly([[yc, xc]], contour):
    """
    # define neighboring points
    xns = xc + np.arange(-radius, radius+1, step=1)
    yns = yc + np.arange(-radius, radius+1, step=1)
    pns = [[yn, xn] for yn in yns for xn in xns] 

    ccontours = []
    for contour in contours:
        if SignedPolygonArea(contour) > radius**2:
            if np.any(skmeasure.points_in_poly(pns, contour)):
                ccontours = ccontours+[contour]
        elif SignedPolygonArea(contour) > areallimit_ctr:
            if np.absolute(contour[:, 0]-yc).min() < radius: 
                if np.absolute(contour[:, 0]-xc).min() < radius: 
                    if np.any(skmeasure.points_in_poly(pns, contour)):
                        ccontours = ccontours+[contour]
    return ccontours


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
    if not ispolys(contours):
        raise ValueError("input contours should be polys")
    if not ispolys(outcontours):
        raise ValueError("input outcontours should be polys")
    incontours = []
    for contour in contours:
        isinsides = np.array([isinside_polygon(contour, outcontour) for outcontour in outcontours])
        isequals = np.array([isequal_polygon(contour, outcontour) for outcontour in outcontours])
        if any(isinsides) and all(~isequals): 
            incontours = incontours+[contour]
    return incontours


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
    return np.all(skmeasure.points_in_poly(poly1,poly2))


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
        area+=SignedPolygonArea(polys[i])
    return area

def SignedPolygonArea(poly):
    """ 
    Calculate the signed area inside a polygon using Shoelace formula. 
    Positive for clock wise (on x-y plane) polygons. 

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
    return area


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
    d=listpoints.dot(unitvec)
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
    thetas=np.linspace(0., 180., num=15, endpoint=False)
    dfs=FeretDiameter(sth,thetas)
    theta_init=thetas[np.argmax(dfs)]
    tomin = lambda theta, sth: -FeretDiameter(sth, theta)

    # fit
    results=minimize(tomin, theta_init,method='Powell', args=(sth,))
    if not results['success']: raise NameError('minimization failed')

    # return
    theta_max=results['x'].sum()
    df_max=FeretDiameter(sth,theta_max)
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



# def FeretDiameter(poly,theta):
#     """
#     Return the Feret's diamter along the direction of theta (y of x.)

#     Parameters
#     ----------
#     poly : (N*2) ndarray of double
#         list of corner coordinates (r,c) of the polygon.
#     theta: (float)
#         the angle of the direction being measured. y of x in degrees. 

#     Returns
#     ----------
#     d_feret: (float)
#         the Feret's diameter in pix
#     """
#     # define unit vector (r,c)=(y,x) along the direction of theta
#     unitvec=np.array([np.sin(np.radians(theta)),np.cos(np.radians(theta))])
#     return np.max(poly.dot(unitvec),axis=0)-np.min(poly.dot(unitvec),axis=0)

# =========== deprecated functions
# def find_isocontours(img, threshold, tocomplete=True, nx=64, ny=64):
#     """
#     Return all isophote contours. Each contours is (N*2) ndarray of double. 
#     High contours are counter-clockwise points and low contours are clockwise 
#     points.

#     if tocomplete is True, for contours that end at edges, add corner points 
#     to complete the loop.
#     """
#     contours = skmeasure.find_contours(img, threshold, fully_connected='low', positive_orientation='high')

#     if tocomplete: # to complete contours
#         for i, contour in enumerate(contours):
#             if SignedPolygonArea(contour)>0: contourtype='high'
#             else: contourtype='low'
#             contourcomplete=complete_contouredges(contour, contourtype=contourtype, nx=nx, ny=ny)
#             contours[i]=contourcomplete

#     return contours



# def complete_contouredges(contour, contourtype='high', nx=64, ny=64):
#     """
#     if the two ends of the contour are on different edges, add a corner(s) to 
#     complete the loop, such that it will still be a closed high/low contour. 
#     """
#     def point_on_edgei(point, nx=64, ny=64):
#         """
#         return which edge is the point on
#         point = (y, x)

#         x = 0 : edgei = 0
#         y = 0 : edgei = 1
#         x = n : edgei = 2
#         y = n : edgei = 3

#         x = 0, y = 0 : edgei = -0.1
#         x = n, y = 0 : edgei = -1.1
#         x = n, y = n : edgei = -2.1
#         x = 0, y = n : edgei = -3.1

#         otherwise np.nan
#         """
#         y, x = point
#         ex=nx-1
#         ey=ny-1

#         if x in [0, ex] or y in [0, ey]:
#             # on edge or corner
#             if x in [0, ex] and y in [0, ey]:
#                 # on corner
#                 if x == 0 and y == 0: edgei = -0.1
#                 elif x == ex and y == 0: edgei = -1.1
#                 elif x == ex and y == ey: edgei = -2.1
#                 elif x == 0 and y == ey: edgei = -3.1
#             else: 
#                 # on edge
#                 if x == 0: edgei = 0
#                 elif y == 0: edgei = 1
#                 elif x == ex: edgei = 2
#                 elif y == ey: edgei = 3            
#         else: 
#             # not on edge nor corner
#             edgei=np.nan
#         return edgei

#     # start complete_highcontouredges() function
#     if contourtype == 'high':
#         p0, p1 = contour[0], contour[-1]
#     elif contourtype == 'low':
#         p0, p1 = contour[-1], contour[0]

#     ei0=point_on_edgei(p0, nx=nx, ny=ny)
#     ei1=point_on_edgei(p1, nx=nx, ny=ny)

#     if ei0 >= 0 and ei1 >= 0: # both on edges
#         if ei0 != ei1: # on different edges
#             # initialization
#             corneriruler=[0, 1, 2, 3, 0, 1]
#             pcorner_list=np.array([[0, 0], [0, nx-1], [ny-1, nx-1], [ny-1, 0]])
#             if ei1 < ei0: ei1 = ei1 + 4
#             corneris=corneriruler[ei0:ei1]
#             pcorners_toadd=pcorner_list[corneris]
#             # adding corners
#             contour=np.concatenate((contour,pcorners_toadd),axis=0)
#             return contour
#         else: 
#             return contour
#     else: 
#         return contour
