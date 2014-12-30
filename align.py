# colordisp/align.py
# ipython -pylab
from pylab import *

import os
import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import ICRS

from scipy.ndimage import interpolation

import tableio

#=================== toolkits

def isInt(val):
    return val == int(val)

#	return tab['f1'][tab['f0']==paramname][0]

def getpixcrd_shift(w1, w2, worldcrd_tar_ref):
	# find the diff of pix coordinates of w1 and w2 wcs systems at poit worldcrd_tar_ref
	pixcrd_shift=w1.wcs_world2pix(worldcrd_tar_ref,0)-w2.wcs_world2pix(worldcrd_tar_ref,0)
	return pixcrd_shift #[col, row]


def cutstampImage(image, xcenter, ycenter, xwidth, ywidth):
	# cut a stamp image from image
	dx=(xwidth-1.)/2.
	dy=(ywidth-1.)/2.
	# check if numbers are integer
	isintegers=[isInt(val) for val in [xcenter,ycenter,dx,dy]]
	if not all(isintegers):
		raise NameError("align.cutstampImage : non-integers received")
	# return stamp image
	return image[xcenter-dx:xcenter+dx+1, ycenter-dy:ycenter+dy+1]


def getalignedImages(obj, bands=('g','r','i'), band_rf='r'):
	filename_rf=obj.dir_obj+'frame-'+band_rf+'.fits'
	w_ref=wcs.WCS(filename_rf)
	worldcrd_obj=np.array([[obj.sdss.ra,obj.sdss.dec]])
	
	# define images_aligned
	header = fits.getheader(filename_rf)
	images_aligned = np.zeros([len(bands),header['NAXIS1'],header['NAXIS2']]) # [band, col, row]
	# read and shift images
	for nb in range(len(bands)):
		filename_sh=obj.dir_obj+'frame-'+bands[nb]+'.fits'
		w_sh=wcs.WCS(filename_sh)
		pixcrd_shift=getpixcrd_shift(w_sh,w_ref,worldcrd_obj)
		image_shi=fits.getdata(filename_sh).transpose()
		images_aligned[nb]=interpolation.shift(image_shi,-pixcrd_shift[0])
	return images_aligned


def getalignedstampImages(obj, bands=('g','r','i'), band_rf='r', xwidth=57, ywidth=57,savefits=True):
	#=== get alinged (w.r.t. reference ('r') band) stamp (xwdith*ywidth) images
	# get center pix coord
	xcenter,ycenter = round(obj.sdss.photoobj['colc']), round(obj.sdss.photoobj['rowc'])
	images_aligned = getalignedImages(obj, bands, band_rf)
	images_aligned_stamp = np.zeros([len(bands),xwidth,ywidth])
	for nb in range(len(bands)):
		# make stamps
		images_aligned_stamp[nb]=cutstampImage(images_aligned[nb], xcenter, ycenter, xwidth, ywidth)#.swapaxes(0,1)
	images_aligned_stamp=np.amax([images_aligned_stamp,np.zeros(images_aligned_stamp.shape)],axis=0)
	# store fits file
	#w_ref=getstampwcs(obj,band_rf,xwidth,ywidth)
	#header = w_ref.to_header()
	if savefits:
		images_aligned_stamp_fits=swapaxes(images_aligned_stamp,1,2)
		for nb in range(len(bands)):
			header=getstampheader(obj,bands[nb],band_rf,xwidth,ywidth)
			filename=obj.dir_obj+'stamp-'+bands[nb]+'.fits'
			prihdu = fits.PrimaryHDU(images_aligned_stamp_fits[nb], header=header)
			prihdu.writeto(filename)
	return images_aligned_stamp #.swapaxes(0,1).swapaxes(1,2)


def getstampwcs(obj, band_rf='r', xwidth=57, ywidth=57):
	# get wcs
	filename_rf=obj.dir_obj+'frame-'+band_rf+'.fits'
	w_ref=wcs.WCS(filename_rf)
	# image width => pix coord of ref pixel
	dx=(xwidth-1.)/2.
	dy=(ywidth-1.)/2.
	w_ref.wcs.crpix = np.array([dx,dy])
	# world coord of reference pixel
	xcenter,ycenter = round(obj.sdss.photoobj['colc']), round(obj.sdss.photoobj['rowc'])
	w_ref.wcs.crval = w_ref.wcs_pix2world([[xcenter,ycenter]], 1)[0]
	return w_ref


def getstampheader(obj, band, band_rf='r', xwidth=57, ywidth=57):
	# get wcs
	filename_rf=obj.dir_obj+'frame-'+band_rf+'.fits'
	#header_rf=fits.getheader(filename_rf)
	filename=obj.dir_obj+'frame-'+band+'.fits'
	header=fits.getheader(filename)
	
	w_ref=wcs.WCS(filename_rf)
	xcenter,ycenter = round(obj.sdss.photoobj['colc']), round(obj.sdss.photoobj['rowc'])
	worldcrds = w_ref.wcs_pix2world([[xcenter,ycenter]], 1)[0]

	header['NAXIS1']=xwidth
	header['NAXIS2']=ywidth
	header['CRPIX1']=(xwidth-1.)/2.
	header['CRPIX2']=(ywidth-1.)/2.
	header['CRVAL1']=worldcrds[0]
	header['CRVAL2']=worldcrds[1]
	
	return header





#prihdu = fits.PrimaryHDU(images_stamp[0], header=header)

#def checkCoords(root, bands, framenames, nf_ref=1):
#	# purpose:
#	# various tests to check the wcs of frame images is the same as sdss
#	# usage:
#	#
#	# root: string, path of the root directory where frame fits files live
#	# framenames: array of strings, paths of the frames align
#	# nf_ref: integer, nth frame is the reference image, default=1
#
#	#==== cd to root directory
#	if bands[nf_ref] != 'r': raise NameError("align.checkCoords: reference frame is not set to r-band")
#	initdir = os.getcwd()
#	if os.path.isdir(root): os.chdir(root)
#	else: raise NameError("align.checkCoords: root directory does not exist")
#	
#	#==== reading ra, dec, colc, rowc from the PhotoObj table
#	worldcrd_obj=array([[getParams('PhotoObj.csv','ra'),getParams('PhotoObj.csv','dec')]]) # in units of deg
#	worldcrdErr_obj=1./3600.*array([[getParams('PhotoObj.csv','raErr'),getParams('PhotoObj.csv','decErr')]]) # in units of deg
#
#	print "worldcrd_obj", worldcrd_obj
#	
#	colcs=[getParams('PhotoObj.csv','colc_'+band) for band in bands]
#	rowcs=[getParams('PhotoObj.csv','rowc_'+band) for band in bands]
#
#	#==== read wcs coordinates from reference frame fits
#	# in wcs_world2pix(X,0) system the corner pixel is (0,0)
#	# in sdss				system the corner pixel is (0.5,0.5)
#	# pixcrd				system mimic sdss system => pix = wcs + 0.5
#	
#	wcspix_sdssorigin =-0.5
#	w_ref=wcs.WCS(framenames[nf_ref])
#	
#	pixcrd_obj_ref=w_ref.wcs_world2pix(worldcrd_obj,0)-wcspix_sdssorigin
#	
#	pixcrdErr_obj_ref=w_ref.wcs_world2pix(worldcrd_obj+worldcrdErr_obj,0)-wcspix_sdssorigin-pixcrd_obj_ref
#	
#	#==== check offset between wcs coordinate from fits and sdss coordinate
#	pixcrd_offset=pixcrd_obj_ref - [[colcs[nf_ref],rowcs[nf_ref]]]
#	print "pixcrd_offset =", pixcrd_offset
#	print "pixcrdErr_obj_ref =", pixcrdErr_obj_ref
#	if any(pixcrd_offset > absolute(pixcrdErr_obj_ref)):
#		raise NameError("align.checkCoords: wcs from frame*.fits is off from sdss coordinate")
#	#==== cd back to initial working directory
#	os.chdir(initdir)


	



