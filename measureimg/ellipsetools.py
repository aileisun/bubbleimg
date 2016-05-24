# ellipsetool.py
# ALS 2016/04/26

"""
Contains tool functions for ellipses parameters and covaraince matrix conversion
"""

import numpy as np

def cov_from_axes(sigma_major, sigma_minor, theta):
    """
    PURPOSE: 
        calcualate the covariance matrix given the axes and orientation 
        of the bivariate Gaussian distribution. 

    PARAMS: 
        sigma_major (length in pixels)
        sigma_minor (length in pixels)
        theta (orientation in degree) y of x

    RETURN:
        cov (array of size 2*2): covariance matrix 

    DESCRIPTION: 
        cov is rotated from [[sigma_major, 0], [0, sigma_minor]]
    """
    theta_rad=np.radians(theta)

    R=np.matrix([[np.cos(theta_rad), -np.sin(theta_rad)], [np.sin(theta_rad), np.cos(theta_rad)]])

    cov=np.matrix([[sigma_major**2,0.],[0.,sigma_minor**2]] )
    cov_rot=R.dot(cov).dot(R.T)

    return cov_rot



def axes_from_cov(cov):
    """
    PURPOSE: 
        calcualate the ellipse parameters: 
            (sigma_major, sigma_minor, theta, eccentricity)
        from the covariance matrix of a bivariate Gaussian distribution. 

    PARAMS: 
        cov (array of size 2*2): covariance matrix 

    RETURN:
        sigma_major (length in pixels)
        sigma_minor (length in pixels)
        theta (orientation in degree) y of x
        # eccen : eccentricity

    DESCRIPTION: 
        The major and minor axis are the eigen values of the cov matrix. And 
        the oreintation theta is from the orientation of the eigen vectors. 
        The eccentricity is defined as sqrt(1- (minor/major)**2)

    # alternative theta formula
    # theta=np.degrees(0.5*np.arctan2(2.*cov[0,1]/(cov[0,0]-cov[1,1])))
    # theta=-np.degrees(np.arctan(vs[1,0]/vs[0,0]))

    """
    # from numpy import linalg

    # eigen decompose the covariance matrix
    [l1,l2],vs=np.linalg.eig(cov)
    l1,l2 = np.abs([l1,l2])    

    # make the larger eigenvalue lambda1 and its vector (major axis) v
    if l1 > l2: 
        lambda1,lambda2=l1,l2
        v=vs[:,0]
    else:
        lambda1,lambda2=l2,l1
        v=vs[:,1]

    # translate to the ellipse proepreties
    sigma_major, sigma_minor = np.sqrt([lambda1,lambda2])
    
    theta=np.degrees(np.arctan2(v[1],v[0]))
    if theta<0.: theta=theta+180.

    # eccen=np.sqrt(1-lambda2/lambda1)

    return sigma_major, sigma_minor, theta #, eccen


