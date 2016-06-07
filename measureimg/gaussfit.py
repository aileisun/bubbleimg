# gaussfit.py
# ALS 2016/04/26

"""
fit bivariate Gaussian distribution to image. 
"""

import numpy as np
import astropy.units as u
from astropy.table import Table
import ellipsetools


def gaussfit_ellipseparams(img, sigma, tofitoffset=False):
    """
    PURPOSE: 
        find best fit 2d Gaussian params from an image
    PARAMS: 
        img (2d np array)
        sigma (float): image noise level
        tofitoffset=False (bool): If True, the fitted image has a freedom of
                                  an overall offset. 
    RETURN:
        m0: sum of pixel values
        xc: x centroid
        yc: y centroid
        a (length in pixels)
        b (length in pixels)
        theta (orientation in degree)
        (offset (float)) optional, will appear if tofitoffset=True. 
    DESCRIPTION: 
        see mom_cov and axes_from_cov
    """

    from scipy.optimize import minimize
    import moment

    nlnlike = lambda *stuffs: -lnlike(*stuffs)

    # use moment params to initialize
    params_init=list(moment.mom_ellipseparams(img))

    if tofitoffset:
        offset_init=0.
        params_init=params_init+[offset_init]

    # fit
    results=minimize(nlnlike, params_init,method='Powell', args=(img, sigma,))
    if not results['success']: raise NameError('gaussfit not successful')
    params_bfit=results['x']

    # getting the order of the axis and PA correct
    [l1,l2]=np.abs(params_bfit[3:5])

    if l1 > l2: 
        a,b=l1,l2
        theta=params_bfit[5]
    else:
        a,b=l2,l1
        theta=params_bfit[5]-90.
    if theta>180.: theta=theta-180.
    params_bfit[3:6]=a,b,theta

    return params_bfit

def img2DGaussian_cov(m0, xc, yc, cov, offset=0., nx=64, ny=64):
    """
    Return an image of bivaraite gaussian from its covaraince

    PARAMS: 
        m0, xc, yc, cov (2*2 np array), offset=0., nx=64, ny=64
    RETURN: img (np array of size nx*ny)
    """
    from scipy.stats import multivariate_normal
    # make an 2d array of [y,x] 
    yx=np.array([[[i,j] for i in range(int(ny))] for j in range(int(nx))]) 
    # make image
    img =offset+m0* multivariate_normal.pdf(yx, mean=[xc,yc], cov=cov)
    return img

def img2DGaussian_axes(m0, xc, yc, a, b, theta, offset=0., nx=64, ny=64):
    """
    Return an image of bivaraite gaussian from its axes

    PARAMS: 
        m0, xc, yc, a, b, theta, offset=0., nx=64, ny=64

    RETURN: img (np array of size nx*ny)
    """
    cov = ellipsetools.cov_from_axes(a, b, theta)
    return img2DGaussian_cov(m0, xc, yc, cov, offset, nx, ny)


def lnlike(params,img_data,sigma):
    """
    PURPOSE: 
        Return the log likelihood of between img_data and img_model, 
        where img_model is the 2D Gaussian specified by params

    PARAMS:
        params (list): (m0, xc, yc, a, b, theta)
        img_data (2d np array)
        sigma (float)
    RETURN: lnlike (float)
    """
    # print type(params),params.shape

    # params=np.append(params) # add in nx, ny to params
    kwargs={'nx':img_data.shape[0],'ny':img_data.shape[1]}
    img_model=img2DGaussian_axes(*params,**kwargs)
    inv_sigma2 = 1.0/(sigma**2)
    lnL=-0.5*np.sum( inv_sigma2*(img_data-img_model)**2 - np.log(inv_sigma2))
    return lnL





# # def sandbox():
#     n=58
#     plt.close()
#     plt.clf()
#     plt.imshow(img,interpolation='nearest',origin='lower left',aspect='equal')
#     # plt.plot([0,n],[0,n],color='black')
#     plt.xlim(0,n)
#     plt.ylim(0,n)
#     plt.show(block=False)
#     n=64
#     plt.close('all')
#     plt.figure()
#     plt.clf()
#     plt.imshow(img_data,interpolation='nearest',origin='lower left',aspect='equal')
#     # plt.plot([0,n],[0,n],color='black')
#     plt.xlim(0,n)
#     plt.ylim(0,n)
#     plt.colorbar()
#     plt.show(block=False)

#     plt.figure()
#     plt.clf()
#     plt.imshow(img_bfit,interpolation='nearest',origin='lower left',aspect='equal')
#     # plt.plot([0,n],[0,n],color='black')
#     plt.xlim(0,n)
#     plt.ylim(0,n)
#     plt.colorbar()
#     plt.show(block=False)
#     plt.figure()
#     plt.clf()
#     plt.imshow(img_data-img_bfit,interpolation='nearest',origin='lower left',aspect='equal')
#     # plt.plot([0,n],[0,n],color='black')
#     plt.xlim(0,n)
#     plt.ylim(0,n)
#     plt.colorbar()
#     plt.show(block=False)

