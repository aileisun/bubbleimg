# psf.py

import numpy as np
from skimage.util import random_noise
# from bubblepy.measureimg import gaussfit
# from scipy.ndimage.filters import convolve


def add_gaussnoise(img, noise_sigma):
    """ add iid gaussian noise on image"""
    if noise_sigma==0:
        return img 
    else:         
        return img+random_noise(np.zeros(img.shape), mode='gaussian',var=noise_sigma**2,clip=False)


# def convolve_gausspsf(img, psf_FWHM, pixelsize=0.396): 
#     """
#     convolve the img with a psf of specified size (sigma) 
    
#     Parameters
#     ------------
#     img: array

#     psf_FWHM
#         [arcsec]

#     pixelsize
#         [arcsec]

#     """
#     psf_sigma=psf_FWHM/2.35482
#     psf_sigma_pix=psf_sigma/pixelsize
#     nx,ny=img.shape
#     psf=psf_gaussian(psf_sigma_pix, nx=nx, ny=ny)
#     result=convolve_psf(img, psf)
#     return result


# def psf_gaussian(sigma, nx=64, ny=64):
#     """
#     return a 2d model of symmetric gaussian PSF, with size specified by sigma
#     sigma is in units of pixel. 
#     """
#     xc=nx*0.5
#     yc=ny*0.5

#     kwargs={'theta':0., 'offset':0., 'nx':nx, 'ny':ny}
#     psf = gaussfit.img2DGaussian_axes(m0=1., xc=xc, yc=yc, a=sigma, b=sigma, **kwargs)
#     return psf


# def convolve_psf(img, psf):
#     """
#     convolve img with psf. Both are 2d array. 
#     notice: psf should centered on    xc=nx*0.5, yc=ny*0.5
#     """
#     result = convolve(img, psf, mode='nearest', origin=0)
#     return result

