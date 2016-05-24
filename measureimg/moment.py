# moments.py
# ALS 2016/04/22

"""
make 2nd moment measurements
"""

import numpy as np
import ellipsetools


def mom_ellipseparams(img):
    """
    PURPOSE: 
        calcualate the ellipse parameters from an image 
        combining the function of two functions: 
            mom_cov
            axes_from_cov
    PARAMS: 
        img
    RETURN:
        m0: sum of pixel values
        xc: x centroid
        yc: y centroid
        a (length in pixels)
        b (length in pixels)
        theta (orientation in degree)

    DESCRIPTION: 
        see mom_cov and axes_from_cov
    """
    from numpy import linalg

    m0, xc, yc, cov = mom_cov(img)
    a, b, theta = ellipsetools.axes_from_cov(cov)

    # return 
    return [m0, xc, yc, a, b, theta]



def mom_cov(img):
    """
    PURPOSE: 
        calcualte the 'key' moments from the image: 
            m0, xc, yc, cov
    PARAMS: 
        img: 2d np array of no units

    RETURN: 
        m0, xc, yc, cov

            m0: sum(I), total flux, sum of pixel values
            xc: sum(x * I)/sum(I), intensity weighted x centroid in pix
            yc: sum(y * I)/sum(I), intensity weighted y centroid in pix
            cov: covariance matrix from 2nd order moments

                 cov(x,x) cov(x,y)      mn[2,0]  mn[1,1]
                 cov(y,x) cov(y,y)  =   mn[1,1]  mn[0,2] 

                 where: mn are central moments normalized by total flux (m0)

    REFERENCES: 
        https://en.wikipedia.org/wiki/Image_moment
        http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/

    """     
    from skimage import measure

    # calculate raw moments, m[ix,iy]
    m = measure.moments(img,order=1)
    # calculate total flux and centroids
    m0=m[0,0] 
    xc = m[1,0]/m0
    yc = m[0,1]/m0  

    #== infer ellips properteis from central 2nd order moments
    # calculate central moments (mu in wiki), notice row = y, column=x.
    mc = measure.moments_central(img,cr=yc,cc=xc,order=2)
    mn=mc/m0 # normalized moments (mu' in wiki)
    cov=np.array([[mn[2,0],mn[1,1]],[mn[1,1],mn[0,2]]])

    return [m0, xc, yc, cov]

