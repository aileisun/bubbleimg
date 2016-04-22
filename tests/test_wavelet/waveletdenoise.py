# waveletdenoise.py
# ALS 2016/04/06

"""
Playing with pywavelet

On-line resources: 
    pywavelet examples: https://pywavelets.readthedocs.org/en/v0.4.0/regression/index.html
    wavelet browser: http://wavelets.pybytes.com


conclusion: use wavedec and waverec to do smoothing + denoising 
            seems to work the best now. 

next step:  1. play with denoising threshold
            2. make some measurements (think about power missed after denoisng)
            3. play with wavelets
"""

from pylab import *

import pywt
import numpy
from astropy.io import fits
from astropy.table import Table
import os

import noiselevel


def main():

    #==== run haar wavelet
    wname='haar'
    wavelet= pywt.Wavelet(wname)

    if not os.path.isdir(wname):
        os.mkdir(wname)

    # run two thresholds 0.5*sigma and 1*sigma
    run_batch(thres_factor=0.5,wavelet=wavelet,directory=wname+'/threshold_05sigma/',suffix='')
    run_batch(thres_factor=1.,wavelet=wavelet,directory=wname+'/threshold_1sigma/',suffix='')


    #==== run haar wavelet
    wname='coif1'
    wavelet= pywt.Wavelet(wname)

    if not os.path.isdir(wname):
        os.mkdir(wname)

    # run two thresholds 0.5*sigma and 1*sigma
    run_batch(thres_factor=0.5,wavelet=wavelet,directory=wname+'/threshold_05sigma/',suffix='')
    run_batch(thres_factor=1.,wavelet=wavelet,directory=wname+'/threshold_1sigma/',suffix='')



def run_batch(thres_factor, wavelet,directory,suffix):

    # setting
    scaling = 1.e15

    if not os.path.isdir(directory):
        os.mkdir(directory)

    # read in data
    filename='stamp-lOIII5008_I.fits'
    hdu=fits.open(filename)
    img=hdu[0].data*scaling # scale up the data

    # set threshold based on sigma 
    thres,sigma=get_fthreshold(img,factor=thres_factor)

    # do the denoisings
    waveletdenoise_2d(img,thres,wavelet=wavelet,directory=directory,suffix=suffix)

    # write summary
    tsummary=Table([[filename],[sigma],[thres],['factor '+str(thres_factor)],[wavelet.name]],names=['filename','sigma','thres_denoise','thres_denoise_type','wavelet'])
    tsummary.write(directory+'params.csv',format='csv')



def waveletdenoise_2d(img,threshold,wavelet=pywt.Wavelet('haar'),directory='./',suffix=''):
    """
    reference: 
    https://blancosilva.wordpress.com/teaching/mathematical-imaging/denoising-wavelet-thresholding/
    """

    #----- decompose
    # nlevel  = int( floor( log2(img.shape[0]) ) ) # 5
    coeffs = pywt.wavedec2( img, wavelet) #, level=nlevel)

    #----- plot original
    plotsave_img(img, directory+'img')
    plot_coeffs(coeffs, directory+'coeff')

    #----- threshold denoising - soft
    operation='_denoised_soft'
    coeffs_new=threshold_wcoeffs(coeffs,threshold,mode='soft')
    img_new = pywt.waverec2(coeffs_new, wavelet)
    plotsave_img(img_new, directory+'img'+operation+suffix)
    plot_coeffs(coeffs_new, directory+'coeff'+operation+suffix)

    # make clipped version
    operation='_denoised_soft_clipthreshold1'
    img_clip=clip_image(img_new, threshold)
    plotsave_img(img_clip, directory+'img'+operation+suffix)

    operation='_denoised_soft_clipthreshold2'
    img_clip=clip_image(img_new, threshold*2.)
    plotsave_img(img_clip, directory+'img'+operation+suffix)

    #----- threshold denoising - hard
    operation='_denoised_hard'
    coeffs_new=threshold_wcoeffs(coeffs,threshold,mode='hard')
    img_new = pywt.waverec2( coeffs_new, wavelet)
    plotsave_img(img_new, directory+'img'+operation+suffix)
    plot_coeffs(coeffs_new, directory+'coeff'+operation+suffix)

    #-----  smooth
    operation='_smooth'
    coeffs_new=smooth_wcoeffs(coeffs,order=1)
    img_new = pywt.waverec2( coeffs_new, wavelet)    
    plotsave_img(img_new, directory+'img'+operation+suffix)
    plot_coeffs(coeffs_new, directory+'coeff'+operation+suffix)

    #----- smooth and soft denoise
    operation='_denoised_soft_smooth'
    coeffs_new=smooth_wcoeffs(coeffs,order=1)
    coeffs_new=threshold_wcoeffs(coeffs_new,threshold,mode='soft')
    img_new = pywt.waverec2( coeffs_new, wavelet)
    plotsave_img(img_new, directory+'img'+operation+suffix)
    plot_coeffs(coeffs_new, directory+'coeff'+operation+suffix)




def threshold_wcoeffs(coeffs,threshold,mode='soft'):
    """
    Thresholding the wavelet coefficients. 

    params: 
        coeffs (list): wavelet coeffs 
                       [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]
        threshold (float): 
        mode='soft'  : 'soft': soft thresholding
                       'hard': hard thresholding
    return 
        coeffs_new
    """
    from copy import deepcopy

    coeffs_new=deepcopy(coeffs)

    # threshold all the detialed coeffs [1:], but not the approximate [1]
    coeffs_new[1:] = map (lambda x: pywt.threshold(x,threshold,mode=mode),
    coeffs[1:])

    return coeffs_new


def smooth_wcoeffs(coeffs,order=1):
    """
    smooth image by setting detailed coefficients zero

    PARAMS: 
        coeffs (list): wavelet coeffs 
                       [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]
        order=1 : smooth to what order
                  setting the detials coefficients 0 up to that order
    RETURN:  
        coeffs_new
    """
    from copy import deepcopy

    coeffs_new=deepcopy(coeffs)

    for j in range(1,order+1):
        # create a list of three zero arrays
        czero=tuple(zeros(coeffs[-j][i].shape) for i in range(3))
        # to replace detailed coeffs with zeros
        coeffs_new[-j]=czero

    return coeffs_new

def clip_image(img,threshold):
    """
    PURPOSE: set image values below threshold to be zero
    PARAMS: 
        img,threshold
    RETURN: 
        img_new
    """
    from copy import deepcopy
    img_new=deepcopy(img)
    img_new[img_new<threshold]=0.
    return img_new

def plotsave_img(img,filename,vmin=0,vmax=30,tosave=True):
    """
    PURPOSE: 
        plot and save image to filename.pdf and filename.npy (np array)
    PARAMS: 
        img,filename,vmin=0,vmax=30,tosave=True    
    """
    clf()
    imshow(img,interpolation='nearest',origin='lower',aspect='equal',cmap='viridis',vmin=vmin,vmax=vmax,extent=[0,img.shape[0],0,img.shape[0]])
    colorbar()
    if tosave:
        savefig(filename+'.pdf')
        save(filename, img)


def plot_coeffs(coeffs,filename,tosave=True):
    """
    DESCRIPTION: make wavelet coefficient plots and store it as filename.pdf

    PARAMS: 
        coeffs: coefficient returned by pywt.wavedec2, in the form
                of [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]. 

        filename: output filename.pdf
    """
    # setup
    clf()

    nlevel=len(coeffs)-1


    na=coeffs[0].shape[0]
    nds=[coeffs[i+1][0].shape[0] for i in range(nlevel)]

    # create zeros array to store coefficients
    ntot=na+sum(nds)
    arr=zeros([ntot,ntot])


    # assign A coeffs 
    arr[0:na, 0:na]=coeffs[0]


    for l in range(nlevel): # for each level of detailed coeffs

        dx=nds[l]
        x0=na+sum(nds[:l])

        # Horizontal
        arr[x0:x0+dx,0:dx]=coeffs[l+1][0]
        # Vertical
        arr[0:dx,x0:x0+dx]=coeffs[l+1][1]
        # Diagonal
        arr[x0:x0+dx,x0:x0+dx]=coeffs[l+1][2]

    # show the coeffs
    imshow(arr,interpolation='nearest',origin='lower left',aspect='equal',cmap='viridis',vmin=-20,vmax=30,extent=[0,ntot,0,ntot])
    colorbar()

    # show the margins
    for l in range(nlevel): # for each level of detailed coeffs

        dx=nds[l]
        x0=na+sum(nds[:l])

        plot([0,x0+dx],[x0,x0],ls='-',color='black')
        plot([x0,x0],[0,x0+dx],ls='-',color='black')

    xlim(0,ntot)
    ylim(0,ntot)

    if tosave: savefig(filename+'.pdf')


def get_fthreshold(img,factor=1.):
    """
    calculating factor threshold = sigma * factor 
    """
    import noiselevel
    # sigma=Table.read('noiselevel.csv',format='csv')['sigma'][0]
    sigma=noiselevel.getnoiselevel(img,ranges=(-30,30),toplot=False)
   
    thres= sigma*factor
    return thres,sigma


def get_uthreshold(img):
    """
    calculating universal threshold = sigma * sqrt(2*log(n))
    where n is the number of points in data ( is that correct ?)
    """
    import noiselevel
    # sigma=Table.read('noiselevel.csv',format='csv')['sigma'][0]
    sigma = noiselevel.getnoiselevel(img,ranges=(-30,30),toplot=False)
   
    thres = sigma*np.sqrt(2*np.log(img.size))
    return thres, sigma
