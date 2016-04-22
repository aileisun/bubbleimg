# denoise.py
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


def wavelet_1d():
    """
    references: 
    http://www.pybytes.com/pywavelets/regression/dwt-idwt.html
    http://www.pybytes.com/pywavelets/regression/multilevel.html
    """
    w=pywt.Wavelet('haar')

    (phi, psi, x) = w.wavefun(level=1)

    x = array([3, 7, 1, 1, -2, 5, 4, 6, 9, 10, 2, 3])


    cA, cD = pywt.dwt(x, 'haar')


    y=pywt.idwt(cA, cD, 'haar')

    threshold=3
    y_d=pywt.idwt(cA, None, 'haar')

    clf()
    plot(x,color='black', label='original')
    plot(y,color='red', label='reconstructed')
    plot(y_d,color='blue', label='de-noised')

    #=== wavedec
    w=pywt.Wavelet('haar')
    coeffs = pywt.wavedec(x, w)
    y=pywt.waverec(coeffs, w)


    #=== wavepacket
    wp = pywt.WaveletPacket(data=x, wavelet='haar', mode='symmetric')


#==== 2d 
def waveletdenoise_2d():
    """
    reference: 
    https://blancosilva.wordpress.com/teaching/mathematical-imaging/denoising-wavelet-thresholding/
    """
    filename='stamp-lOIII5008_I.fits'
    hdu=fits.open(filename)

    img=hdu[0].data

    img=img*1.e15

    #----- plot original
    fileout='img'
    clf()
    imshow(img,interpolation='nearest',origin='lower')
    savefig(fileout+'.pdf')
    save(fileout, img_new)

    #----- decompose
    wavelet = pywt.Wavelet('haar')
    levels  = int( floor( log2(img.shape[0]) ) ) # 5

    coeffs = pywt.wavedec2( img, wavelet, level=levels)


    #----- threshold denoising - soft
    fileout='img_denoised_soft_1'
    threshold = 1.
    coeffs_new = map (lambda x: pywt.threshold(x,threshold,mode='soft'),
    coeffs)
    img_new = pywt.waverec2( coeffs_new, wavelet)

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig(fileout+'.pdf')
    save(fileout, img_new)


    #----- threshold denoising - hard
    fileout='img_denoised_hard_1'
    threshold = 1.
    coeffs_new = map (lambda x: pywt.threshold(x,threshold,mode='hard'),
    coeffs)
    img_new = pywt.waverec2( coeffs_new, wavelet)

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig(fileout+'.pdf')
    save(fileout, img_new)

    #-----  smooth

    fileout='img_denoised_soft_1_smooth'

    from copy import deepcopy
    coeffs_new=deepcopy(coeffs)

    # make 3 zero arrays
    shape_detailed=coeffs[-1][0].shape
    cHVD_zero=tuple(zeros(shape_detailed) for i in range(3))

    # replace the 1st detailed coeffs with zeros
    coeffs_new[-1]=cHVD_zero

    img_new = pywt.waverec2( coeffs_new, wavelet)

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig('img_smooth.pdf')

    # smooth and soft denoise
    threshold = 1.
    coeffs_new = map (lambda x: pywt.threshold(x,threshold,mode='soft'),
    coeffs_new)

    img_new = pywt.waverec2( coeffs_new, wavelet)

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig(fileout+'.pdf')
    save(fileout, img_new)


def waveletpacket2d():
    """

    reference: 
    http://www.pybytes.com/pywavelets/regression/wp2d.html

    try using waveletpacket instead of wavedec

    experience learned: 
        Whenever WaveletPacket2D is called or when coeffs (.data) of 
        a node is changed, the decomposition is done to update all 
        the down stream levels. 

        e.g., when wp['v'].data is updated, wp['va'], wp['vv'], wp['vh'], wp['vd'], wp['vaa']..., etc. are all updated. 

        When one wants to update the upstream from new downstream coeffs,
        one needs to call 

            node.reconstruct(update=True). 

        At this point, only one level is reconstruction is done, i.e. 

            node.data is updated from the current
                node['a'] node['v'] node['h'] node['d']. 

        At this point, if any of the downstream nodes is not updated 
        with respect to its downstream ndoes, e.g., if node['av'] has been 
        changed but node['a'] has not been updated, then the code will refuse
        to update node before node['a'] is updated. 

    """


    filename='stamp-lOIII5008_I.fits'
    fileout='img_smooth_wpacket'
    hdu=fits.open(filename)

    img=hdu[0].data

    img=img*1.e15

    #----- plot original
    clf()
    imshow(img,interpolation='nearest',origin='lower')
    # savefig('img.pdf')

    # decomposition
    wp=pywt.WaveletPacket2D(img,wavelet='haar',maxlevel=3)

    # kill the first level details
    for node in ['d','h','v']:
        wp[node].data=zeros(wp[node].data.shape)

    # reconstruct
    wp.reconstruct(update=True)
    img_new=wp.data

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig(fileout+'.pdf')

    save(fileout, img_new)

def waveletpacket2d_threshold():
    """
    packet is not good for denoising, it's too flexible in terms of 
    decomposing ....
    """

    #----- setups

    filename='stamp-lOIII5008_I.fits'
    hdu=fits.open(filename)

    img=hdu[0].data
    img=img*1.e15

    threshold=1.

    maxlevel  = int( floor( log2(img.shape[0]) ) ) # 5


    # decomposition
    wp=pywt.WaveletPacket2D(img,wavelet='haar',maxlevel=maxlevel, mode='sym')
    wp['aaaaa']
    # thresholding
    leaf_nodes=[n.path for n in wp.get_leaf_nodes()]
    print leaf_nodes

    for n in leaf_nodes:
        wp[n].data=pywt.threshold(wp[n].data,threshold,mode='soft')

    parent_nodes=unique([wp[n].parent.path for n in leaf_nodes])[::-1]

    for n in ['aaaaa']+list(parent_nodes):
        print n
        wp[n].reconstruct(update=True)
    # .......... got stuck here ..........

    for n in ['aaaaa','aaaa']:
        print n
        wp[n].reconstruct(update=True)


    # reconstruct
    wp.reconstruct(update=True)
    img_new=wp.data

    clf()
    imshow(img_new,interpolation='nearest',origin='lower')
    savefig('img_denoised_soft_1_wpacket.pdf')




