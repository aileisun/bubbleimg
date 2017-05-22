# sdssdisp/alignstamp.py
# ipython -pylab
# ALS reorganized 2015/08/10

"""
PURPOSE: make aligned stamp multi-band images of an SDSS object. 

"""

import os
import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.ndimage import interpolation
from astroquery.sdss import SDSS


def write_alignedstampImages(obj, bands=('g','r','i'), band_rf='r', xwidth=64, ywidth=64, clipnegative=False, overwrite=True):
	"""
	PURPOSE: write aligned stamp images of an SDSS object

	INPUT: 
		obj (obj of class obsobj): an object containing information like RA, Dec, and SDSS PhotoObj. 
			Used attributes:  obj.dir_obj, obj.sdss.ra, obj.sdss.dec, obj.sdss.photoobj
		bands=('g','r','i'): list of bands to use
		band_rf='r':         the reference band. All the other band have images shifted (interpolated) to align with the reference band. 
		xwidth=64 (int):	 width of the image
		ywidth=64 (int):	 height of the image
		clipnegative=False (bool): whether to set negative values to zero

	READ INPUT: 
		obj.dir_obj+'frame-'+bands[nb]+'.fits'

	WRITE OUTPUT: 
		obj.dir_obj+'stamp-'+bands[nb]+'.fits'

	DESCRIPTION: 
		The header of the stamp fits file is taken from the frame fits files, with updated wcs of the reference frame. 
	"""
	#=== get alinged (w.r.t. reference ('r') band) stamp (xwdith*ywidth) images
	# get center pix coord
	print "[alignstamp] aligning and stamping images "+str(bands)
	xcenter,ycenter = round(obj.sdss.photoobj['colc']), round(obj.sdss.photoobj['rowc'])
	images_aligned = getalignedImages(obj, bands, band_rf)
	images_aligned_stamp = np.zeros([len(bands),xwidth,ywidth])
	for nb in range(len(bands)):
		# make stamps
		images_aligned_stamp[nb]=cutstampImage(images_aligned[nb], xcenter, ycenter, xwidth, ywidth)#.swapaxes(0,1)
	if clipnegative:
		images_aligned_stamp=np.amax([images_aligned_stamp,np.zeros(images_aligned_stamp.shape)],axis=0)
	# store fits file
	images_aligned_stamp_fits=np.swapaxes(images_aligned_stamp,1,2)
	for nb in range(len(bands)):
		header=getstampheader(obj,bands[nb],band_rf,xwidth,ywidth)
		filename=obj.dir_obj+'stamp-'+bands[nb]+'.fits'
		prihdu = fits.PrimaryHDU(images_aligned_stamp_fits[nb], header=header)
		if os.path.isfile(filename):
			os.remove(filename)
			prihdu.writeto(filename, overwrite=overwrite)
		else:
			prihdu.writeto(filename, overwrite=overwrite)



def getalignedImages(obj, bands=('g','r','i'), band_rf='r'):
	"""
	PURPOSE: Return n-band aligned frame images of an SDSS object. 
			 This function is used by write_alignedstampImages() as an intermediate 
			 step to get aligned stamp images. 

	NOTE: 
			SDSS Frame Images of the desired bands have to be save to
				obj.dir_obj+'frame-'+bands+'.fits'
			for this funciton to access. 

	INPUT: 
			obj (obj of class obsobj): an object containing information like RA, Dec, and SDSS PhotoObj. 
				of an astro-object with info like RA, Dec, and SDSS PhotoObj. 
				Used attributes: 
					obj.dir_obj, obj.sdss.ra, obj.sdss.dec
			bands=('g','r','i'): list of bands to use
			band_rf='r':         the reference band. All the other band have images shifted (interpolated) to align with the reference band. 

	READ INPUT: 
			Frame images will be read from files: 
			obj.dir_obj+'frame-'+bands+'.fits'


	RETURN OUTPUT: 
			images_aligned (array of shape [nb,ncol,nrow]): aligned frame images of the bands. 
	"""
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



def cutstampImage(image, xcenter, ycenter, xwidth, ywidth, fill_value=0.):
	"""
	Cut out a stamp image from image. If the stamp goes outside the image the edge is filled with fill_value. 

	Parameters
	------
			image (2d array): image to be cut out from 
			xcenter (int): 	  x center of the stamp 
			ycenter (int):    y center of the stamp 
			xwidth  (int):    x width of the stamp
			ywidth  (int):    y width of the stamp
			fill_value (float): value filing the edge if stamp goes outside the image

	Return
	------
	stampimage: 2d array

	Note
	------
	If xwidth, ywidth are even numbers, the center is on [xwidth/2., ywidth/2.]
	of the stamp image. Otherwise if odd, the center is on [(xwidth-1)/2.,(ywidth-1)/2.]
	"""
	dx = np.absolute((xwidth)/2.)
	dy = np.absolute((ywidth)/2.)

	xwidth_img, ywidth_img = image.shape

	if (xwidth % 2 == 0) and (ywidth % 2 == 0):
		epsilon = 0.
	elif (xwidth % 2 == 1) and (ywidth % 2 == 1):
		epsilon = 0.5
	else: 
		raise ValueError("alignstamp.cutstampImage: xwdith and ywidth have different parity")

	# set stamp cut boundary in image coordinate
	x0 = xcenter-dx+epsilon
	x1 = xcenter+dx+epsilon
	y0 = ycenter-dy+epsilon
	y1 = ycenter+dy+epsilon

	# determine the pixel range to be passed over
	x_min_stp = max(0, -x0)
	x_max_stp = min(xwidth, xwidth - (x1 - image.shape[0]))
	y_min_stp = max(0, -y0)
	y_max_stp = min(ywidth, ywidth - (y1 - image.shape[1]))

	x_min_img = max(0, x0)
	x_max_img = min(xwidth_img, x1)
	y_min_img = max(0, y0)
	y_max_img = min(ywidth_img, y1)


	# # check if numbers are integer
	isintegers=[isInt(val) for val in [x0,x1,y0,y1]]
	if not all(isintegers):
		raise ValueError("alignstamp.cutstampImage: non-integers received")

	# stamping
	stampimage = np.ones([xwidth, ywidth])*fill_value

	stampimage[int(x_min_stp):int(x_max_stp), int(y_min_stp):int(y_max_stp)] = image[int(x_min_img):int(x_max_img), int(y_min_img):int(y_max_img)]

	# return stamp image
	return stampimage



#===== toolkits

def isInt(val):
	""" Return (bool) whether the input is a integer """
	return val == int(val)

def getpixcrd_shift(w1, w2, worldcrd_tar_ref):
	"""
	PURPOSE: Find the diff of pix coordinates of w1 and w2 wcs systems at the RA, Dec poit of worldcrd_tar_ref

	INPUT: 
			w1	(object of class WCS): first wsc system 
			w2	(object of class WCS): second wsc system 
			worldcrd_tar_ref (array) : world coordintae point array([[ra, dec]]) 

	OUTPUT: pixcrd_shift (list): [delta column, delta row]
	"""
	pixcrd_shift=w1.wcs_world2pix(worldcrd_tar_ref,0)-w2.wcs_world2pix(worldcrd_tar_ref,0)
	return pixcrd_shift #[col, row]


def getstampheader(obj, band, band_rf='r', xwidth=64, ywidth=64):
	"""
	PURPOSE: Return the header of the aligned stamp image of a specified band. 

	DESCRIPTION: The header template is taken from the band frame file
						obj.dir_obj+'frame-'+band+'.fits'
				 with wcs coordinate center replaced by the object location in the 
				 reference band, and array size updated to the new size. 

	INPUT: 
			obj (object of the class obsobj)
				required atributes: obj.dir_obj, obj.sdss.photoobj
			band 		(string)
			band_rf='r' (string)
			xwidth=64 	(int)
			ywidth=64	(int)

	READ INPUT: 
			obj.dir_obj+'frame-'+band_rf+'.fits'
			obj.dir_obj+'frame-'+band+'.fits'

	OUTPUT: header
	"""
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

