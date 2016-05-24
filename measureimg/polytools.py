# polytools.py
# ALS 2016/04/29

"""
Contain functions for shape measurements on polygons. 

References: 
https://en.wikipedia.org/wiki/Shape_factor_(image_analysis_and_microscopy)
https://en.wikipedia.org/wiki/Feret_diameter
PEER REVIEW Particle Shape Factors and Their Use in Image Analysis Part 1: Theory by Eric Olson
"""

import numpy as np



def NetPolygonArea(polys):
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
    Positive for counter-clock wise polygons. 

    Parameters
    ----------
    poly : (N*2) ndarray of double
        list of corner coordinates (r,c) of the polygon 
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


def FeretDiameter(poly,theta):
    """
    Return the Feret's diamter along the direction of theta (y of x.)

    Parameters
    ----------
    poly : (N*2) ndarray of double
        list of corner coordinates (r,c) of the polygon.
    theta: (float)
        the angle of the direction being measured. y of x in degrees. 

    Returns
    ----------
    d_feret: (float)
        the Feret's diameter in pix
    """
    # define unit vector (r,c)=(y,x) along the direction of theta
    unitvec=np.array([np.sin(np.radians(theta)),np.cos(np.radians(theta))])
    return np.max(poly.dot(unitvec),axis=0)-np.min(poly.dot(unitvec),axis=0)

def FeretD_max(poly):
    """
    Return the maximum Feret's diameter and its orientaion (y of x). 

    Parameters
    ----------
    poly : (N*2) ndarray of double
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
    thetas=np.linspace(0., 180., num=120, endpoint=False)
    dfs=FeretDiameter(poly,thetas)
    theta_init=thetas[np.argmax(dfs)]
    tomin = lambda theta, poly: -FeretDiameter(poly, theta)

    # fit
    results=minimize(tomin, theta_init,method='Powell', args=(poly,))
    if not results['success']: raise NameError('minimization failed')

    # return
    theta_max=results['x'].sum()
    df_max=FeretDiameter(poly,theta_max)
    return df_max, theta_max


def FeretD_min(poly):
    """
    Return the minimum Feret's diameter and its orientaion (y of x). 

    Parameters
    ----------
    poly : (N*2) ndarray of double
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
    thetas=np.linspace(0., 180., num=120, endpoint=False)
    dfs=FeretDiameter(poly,thetas)
    theta_init=thetas[np.argmin(dfs)]
    tomin = lambda theta, poly: FeretDiameter(poly, theta)

    # fit
    results=minimize(tomin, theta_init,method='Powell', args=(poly,))
    if not results['success']: raise NameError('minimization failed')

    # return
    theta_min=results['x'].sum()
    df_min=FeretDiameter(poly,theta_min)
    return df_min, theta_min


def FeretD_max90(poly):    
    """
    Return Feret's diameter along the direction perpendicular to the maximum
    Feret's diameter. 
    """
    df_max, theta_max=FeretD_max(poly)
    theta_max90=theta_max-90.
    df_max90=FeretDiameter(poly,theta_max90)
    return df_max90, theta_max90

def FeretAspectRatio(poly):
    """
    Return the aspect ratio as defined by FeretD_max90/FeretD_max
    where FeretD_max is the maximum Feret's Diamter and FeretD_max90 is the 
    Feret's diameter perpendicular to the maximum's direction. 
    """
    return FeretD_max90(poly)[0]/FeretD_max(poly)[0]


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
    poly=deepcopy(poly)

    xmin,xmax=np.min(poly[:,1],axis=0.),np.max(poly[:,1],axis=0.)
    ymin,ymax=np.min(poly[:,0],axis=0.),np.max(poly[:,0],axis=0.)

    poly[:,1]=poly[:,1]-xmin
    poly[:,0]=poly[:,0]-ymin

    img=np.zeros([int(ymax-ymin),int(xmax-xmin)])

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

# def Elongation(poly):
