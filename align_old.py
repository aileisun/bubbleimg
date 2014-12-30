# colordisp/align.py
# ipython -pylab

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

def plotandsave_test(image,filename='test'):
	clf()
	imshow(log(clip(image.transpose(),0,amax(image))),origin='lower', interpolation='nearest')
	savefig(filename+'.png')

def getParams(filename, paramname):
	# reading params from PhotoObj table
	tab=tableio.rcsv(filename)
	return tab['f1'][tab['f0']==paramname][0]

def getpixcrd_shift(w1, w2, worldcrd_tar_ref):
	# find the diff of pix coordinates of w1 and w2 wcs systems at poit worldcrd_tar_ref
	pixcrd_shift=w1.wcs_world2pix(worldcrd_tar_ref,0)-w2.wcs_world2pix(worldcrd_tar_ref,0)
	return pixcrd_shift #[col, row]

def getworldcrd_obj(root='.'):
	# get ra and dec of the object from PhotoObj table
	path=os.path.join(root,'PhotoObj.csv')
	ra=getParams(path,'ra')
	dec=getParams(path,'dec')
	return np.array([[ra,dec]]) # [ra, dec] in units of deg

def getpixcrd_obj(root='.'):
	# get row and col of the object from PhotoObj table
	path=os.path.join(root,'PhotoObj.csv')
	colc=getParams(path,'colc')
	rowc=getParams(path,'rowc')
	return np.array([[colc,rowc]]) # [colc, rowc]

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



#======================= main operations

def getpixShifts(root, bands, framenames, nf_ref=1):
	# get shifts between to-shift (_shi) frame and reference frame (_ref)
	w_ref=wcs.WCS(os.path.join(root,framenames[nf_ref]))
	worldcrd_obj=getworldcrd_obj(root)
	pixcrd_shifts = [getpixcrd_shift(wcs.WCS(os.path.join(root,framename)),w_ref,worldcrd_obj)[0] for framename in framenames]
	if any(pixcrd_shifts[nf_ref] != array([0.,0.])) :
		raise NameError("align.getpixShifts: pixcrd_shifts[nf_ref] =! 0")
	return getpixcrd_shifts #[nf, col, row]

def getalignedImages(root, bands, framenames, nf_ref=1):
	# get images of frames aligned w.r.t. the reference ('r') frame
	w_ref=wcs.WCS(os.path.join(root,framenames[nf_ref]))
	worldcrd_obj=getworldcrd_obj(root)

	# define images_aligned
	header = fits.getheader(os.path.join(root,framenames[nf_ref]))
	images_aligned = np.zeros([len(bands),header['NAXIS1'],header['NAXIS2']]) # [band, col, row]
	# read and shift images
	for nf in range(len(bands)):
		w_shi=wcs.WCS(os.path.join(root,framenames[nf]))
		pixcrd_shift=getpixcrd_shift(w_shi,w_ref,worldcrd_obj)
		image_shi=fits.getdata(os.path.join(root,framenames[nf])).transpose()
		images_aligned[nf]=interpolation.shift(image_shi,-pixcrd_shift[0])
	return images_aligned

def getalignedstampImages(root, bands, framenames, nf_ref=1, xwidth=57, ywidth=57):
	# get alinged (w.r.t. reference ('r') band) stamp (xwdith*ywidth) images
	xcenter,ycenter = np.around(getpixcrd_obj(root)[0])
	images_aligned = getalignedImages(root, bands, framenames, nf_ref=1)
	images_aligned_stamp = np.zeros([len(bands),xwidth,ywidth])
	for nf in range(len(bands)):
		images_aligned_stamp[nf]=cutstampImage(images_aligned[nf], xcenter, ycenter, xwidth, ywidth)
	images_aligned_stamp=np.amax([images_aligned_stamp,np.zeros(images_aligned_stamp.shape)],axis=0)
	return images_aligned_stamp #.swapaxes(0,1).swapaxes(1,2)



def checkCoords(root, bands, framenames, nf_ref=1):
	# purpose:
	# various tests to check the wcs of frame images is the same as sdss
	# usage:
	#
	# root: string, path of the root directory where frame fits files live
	# framenames: array of strings, paths of the frames align
	# nf_ref: integer, nth frame is the reference image, default=1

	#==== cd to root directory
	if bands[nf_ref] != 'r': raise NameError("align.checkCoords: reference frame is not set to r-band")
	initdir = os.getcwd()
	if os.path.isdir(root): os.chdir(root)
	else: raise NameError("align.checkCoords: root directory does not exist")
	
	#==== reading ra, dec, colc, rowc from the PhotoObj table
	worldcrd_obj=array([[getParams('PhotoObj.csv','ra'),getParams('PhotoObj.csv','dec')]]) # in units of deg
	worldcrdErr_obj=1./3600.*array([[getParams('PhotoObj.csv','raErr'),getParams('PhotoObj.csv','decErr')]]) # in units of deg

	print "worldcrd_obj", worldcrd_obj
	
	colcs=[getParams('PhotoObj.csv','colc_'+band) for band in bands]
	rowcs=[getParams('PhotoObj.csv','rowc_'+band) for band in bands]

	#==== read wcs coordinates from reference frame fits
	# in wcs_world2pix(X,0) system the corner pixel is (0,0)
	# in sdss				system the corner pixel is (0.5,0.5)
	# pixcrd				system mimic sdss system => pix = wcs + 0.5
	
	wcspix_sdssorigin =-0.5
	w_ref=wcs.WCS(framenames[nf_ref])
	
	pixcrd_obj_ref=w_ref.wcs_world2pix(worldcrd_obj,0)-wcspix_sdssorigin
	
	pixcrdErr_obj_ref=w_ref.wcs_world2pix(worldcrd_obj+worldcrdErr_obj,0)-wcspix_sdssorigin-pixcrd_obj_ref
	
	#==== check offset between wcs coordinate from fits and sdss coordinate
	pixcrd_offset=pixcrd_obj_ref - [[colcs[nf_ref],rowcs[nf_ref]]]
	print "pixcrd_offset =", pixcrd_offset
	print "pixcrdErr_obj_ref =", pixcrdErr_obj_ref
	if any(pixcrd_offset > absolute(pixcrdErr_obj_ref)):
		raise NameError("align.checkCoords: wcs from frame*.fits is off from sdss coordinate")
	#==== cd back to initial working directory
	os.chdir(initdir)


	




#	for nf in range(len(bands)):
#		print "band", bands[nf], "has shift (col, row): ", pixcrd_shifts[i]
#	if any(pixcrd_shifts[nf_ref] != array([0.,0.])) :
#		raise NameError("align: pixcrd_shifts[nf_ref] =! 0")
#
#
#	for nf in range(len(bands)):
#
#
#	
#	image_ref=fits.getdata(framenames[nf_ref])
#	plotandsave_test(image_ref,'test_ref')
#
#	for nf in range(len(bands)) :
#		w_shi=wcs.WCS(framenames[nf])
#
#		pixcrd_shift=getpixcrd_shift(w_shi,w_ref,worldcrd_obj)
#
#		image_shi=fits.getdata(framenames[nf])
#		plotandsave_test(image_shi,'test_shi')
#
#		image_shifted=interpolation.shift(image_shi,-roll(pixcrd_shift[0],1))
#		plotandsave_test(image_shifted,'test_shifted')
#
#

	

	
	
	
#	
#	# temporary using the first framenames as example
##	framename=framenames[0]
#
#	framefits=fits.open(os.path.join(root,framename))
#	
##	framefits.info()
#
#	w=wcs.WCS(framefits[0].header)
#
#	pixcrd = np.array([[50,50],[100,300],[2000,1000]])
#
#	worldcrd=w.wcs_pix2world(pixcrd,1)
#
#
#	c = ICRS(ra=worldcrd[0,0], dec=worldcrd[0,1],unit=(u.degree, u.degree))
#
#	c.to_string()
#
#
#	pixcrd2=w.wcs_world2pix(worldcrd,1)
#
##framename='/Users/aisun/Documents/Astro/Thesis/bbselection/data/SDSS/J2240-0927/frame-r-002576-3-0124.fits'
#root='/Users/aisun/Documents/Astro/Thesis/bbselection/data/SDSS/J2240-0927/'
