# main.py
# ALS 2016/05/03

"""
given obsobj, make all the images

"""

import alignstamp
import imagedisp_util
import makemap

def obj_makeblobmaps(obj, bandline='r',bandconti='z',update=False):
    """
    For a given obsobj object, make all the images, including: 

    alignstamp:
        frame-[g-z].fits
        stamp-[g-z].fits
        images_stamp.npy

    imagedisp_util:
        HumVI_gri.png

    makemap:
        spec.fits
        spec.pdf
        stamp-r-uconti.fits, stamp-r-gconti.fits, stamp-r-iconti.fits, stamp-r-zconti.fits 
        stamp-lOIII5008_F.fits
        stamp-lOIII5008_I.fits
        stamp-lOIII5008_I.pdf
        stamp-lOIII5008_I.png


    Paramters
    ------
     obj

    Write Output
    ------
    obj.dir_obj/...
    """

    # make stamps if they do not exist already
    alignstamp.objw_makeall(obj, update=update)

    # make color images
    imagedisp_util.objw_HumVIgriimages(obj, update=update)
    
    # subtract stamps and make OIII maps
    makemap.objw_makeall(obj, bandline=bandline,bandconti=bandconti,update=update)
