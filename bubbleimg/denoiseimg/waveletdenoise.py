# waveletdenoise.py
# ALS 2016/05/03

"""
denoise images using wavelet
"""

import numpy as np
import matplotlib.pyplot as plt
import pywt


def waveletdenoise_2d(img,threshold,wavelet='coif1',mode='soft',level=3,smooth=False):
    """
    Denoise an image using wavelet

    Parameters
    ------
    img: 2d np array
    threshold: float
    wavelet='coif1' : string
        can be 'haar', 'coif1', etc. see pywt documentations. 
    mode='soft': string
        can be 'soft', 'hard', 'none' (for no threshold denoising)
    level=3: int
        number of levels to do wavelet decomposition. Thresholding will be 
        applied to all levels. 
    smooth=False: bool

    Reference
    ------
    https://blancosilva.wordpress.com/teaching/mathematical-imaging/denoising-wavelet-thresholding/

    Description
    ------
    By defailt, wavelet decomposition is done to the coarest level and 
    thresholding is applied uniformly to all levels of detailed coefficients. 

    If smooth == True, then set the finest detailed coeffs zero. 

    Need to modify the code to enable plotting coeffs. 
    """

    w=pywt.Wavelet(wavelet)
    #----- decompose
    # nlevel  = int( floor( log2(img.shape[0]) ) ) # 5
    coeffs = pywt.wavedec2(img, w, level=level)

    #----- threshold denoising - soft
    coeffs_new=threshold_wcoeffs(coeffs,threshold,mode=mode)
    if smooth:
        coeffs_new=smooth_wcoeffs(coeffs_new,order=1)

    img_new = pywt.waverec2(coeffs_new, wavelet=w)

    # optional plotting:
    # plotsave_img(img_new, directory+'img'+operation+suffix)
    # plot_coeffs(coeffs_new, directory+'coeff'+operation+suffix)

    return img_new



#====== tool sets: 

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

    if mode == 'none': 
        return coeffs_new
    else:
        # threshold all the detialed coeffs [1:], but not the approximate [1]
        coeffs_new[1:] = map (lambda x: pywt.threshold(x,threshold,mode=mode),
        coeffs[1:])

        return coeffs_new


def smooth_wcoeffs(coeffs,order=1):
    """
    smooth image by setting the finest detailed coefficients zero

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
        czero=tuple(np.zeros(coeffs[-j][i].shape) for i in range(3))
        # to replace detailed coeffs with zeros
        coeffs_new[-j]=czero

    return coeffs_new



#======== plotting tools ...... need to add modulo names

def plotsave_img(img,filename,vmin=None,vmax=None,tosave=True):
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


def plot_coeffs(coeffs,filename,vmin=None,vmax=None,tosave=True):
    """
    DESCRIPTION: make wavelet coefficient plots and store it as filename.pdf

    PARAMS: 
        coeffs: coefficient returned by pywt.wavedec2, in the form
                of [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]. 

        filename: output filename.pdf
    """
    # setup
    nlevel=len(coeffs)-1
    na=coeffs[0].shape[0]
    nds=[coeffs[i+1][0].shape[0] for i in range(nlevel)]

    # create zeros array to store coefficients
    ntot=na+sum(nds)
    arr=np.zeros([ntot,ntot])
    # assign A coeffs 
    arr[0:na, 0:na]=coeffs[0]
    # assign detailed coeffs
    for l in range(nlevel): # for each level of detailed coeffs
        dx=nds[l]
        x0=na+sum(nds[:l])
        # Horizontal
        arr[x0:x0+dx,0:dx]=coeffs[l+1][0]
        # Vertical
        arr[0:dx,x0:x0+dx]=coeffs[l+1][1]
        # Diagonal
        arr[x0:x0+dx,x0:x0+dx]=coeffs[l+1][2]

    # plotting - coeffs
    clf()
    imshow(arr,interpolation='nearest',origin='lower left',aspect='equal',cmap='viridis',vmin=None,vmax=None,extent=[0,ntot,0,ntot])
    colorbar()

    # plotting - margins
    for l in range(nlevel): # for each level of detailed coeffs
        dx=nds[l]
        x0=na+sum(nds[:l])
        plot([0,x0+dx],[x0,x0],ls='-',color='black')
        plot([x0,x0],[0,x0+dx],ls='-',color='black')

    xlim(0,ntot)
    ylim(0,ntot)

    if tosave: savefig(filename+'.pdf')


# def get_fthreshold(img,factor=1.):
#     """
#     not used right now...
#     calculating factor threshold = sigma * factor 
#     """
#     import noiselevel
#     # sigma=Table.read('noiselevel.csv',format='csv')['sigma'][0]
#     sigma=noiselevel.getnoiselevel(img,ranges=(-30,30),toplot=False)
   
#     thres= sigma*factor
#     return thres, sigma


# def get_uthreshold(img):
#     """
#     not used right now...
#     calculating universal threshold = sigma * sqrt(2*log(n))
#     where n is the number of points in data ( is that correct ?)
#     """
#     import noiselevel
#     # sigma=Table.read('noiselevel.csv',format='csv')['sigma'][0]
#     sigma = noiselevel.getnoiselevel(img,ranges=(-30,30),toplot=False)
   
#     thres = sigma*np.sqrt(2*np.log(img.size))
#     return thres, sigma
