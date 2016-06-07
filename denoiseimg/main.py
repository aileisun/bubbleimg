# main.py
# ALS 2016/05/03

"""
given obsobj, make the denoised images
"""
import os
from astropy.io import fits

import waveletdenoise
import noiselevel

def obj_makedenoised_fits(obj,**kwargs):
    """
    wraper of makedenoised for obj. 
    """
    makedenoised(dir_obj=obj.dir_obj,**kwargs)    

def dir_makedenoised_fits(dir_obj,filename='stamp-lOIII5008_I.fits',wavelet='coif1',mode='soft',smooth=False,thresholdNR=1.,update=False):
    """
    For a given obsobj object and image filename, make denoised image.

    Parameters
    ------
    dir_obj
    filename='stamp-lOIII5008_I.fits'
    wavelet='coif1'
    mode='soft'
    smooth=False
    thresholdNR=1
        threshold to noise ratio, where noise is obtained from gaussian 
        fitting of image pix value histogram. 

    Write Outputs
    ------
    dir_obj/filename+'_denoised.fits'

    """
    # set up
    filein=dir_obj+filename
    fileout=dir_obj+os.path.splitext(filename)[0]+'_denoised.fits'

    if not os.path.isfile(fileout) or update:
        # read in 
        img=fits.getdata(filein)
        header=fits.getheader(filein)

        # operation
        # get noiselevel and write in noiselevel.csv
        sigma=noiselevel.load_noiselevel(dir_obj, filename,update=False)
        threshold = thresholdNR * sigma
        # denoise
        img_new=waveletdenoise.waveletdenoise_2d(img,threshold,wavelet=wavelet,mode=mode,smooth=smooth)

        # write header
        header['HISTORY']='Wavelet denoised by ALS, see WAVELET, WDENOISE, SMOOTHED'
        header['WAVELET']=wavelet
        header['WDENOISE']=mode
        header['SMOOTHED']=str(smooth)

        # write output
        prihdu = fits.PrimaryHDU(img_new, header=header)
        prihdu.writeto(fileout, clobber=True)
    else:
        print "skip dir_makedenoised_fits() as file exist"
