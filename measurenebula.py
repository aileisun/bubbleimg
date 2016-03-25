# measurenebula.py
# ALS 2015/08/14

"""
PURPOSE: measure the nebula properties, like size, area, luminosities from SDSS images, and store them in the measure....txt file under each batch directory.  

"""
from pylab import *
import os
# import numpy as np
import astropy.units as u
from astropy.table import Table,hstack,vstack
from astropy.io import fits

# import class_obsobj




def obj_measureISO_lOIII5008(obj,isophotocut_base=1.e-15*u.Unit('erg s-1 cm-2 arcsec-2'),smoothing=4.,towriteimgmos=True):
	"""
	PURPOSE: measure isophotal quantities of OIII for the obj, and write mosaiced image
	"""
	filein=obj.dir_obj+'stamp-lOIII5008_I.fits'

	if not os.path.isfile(filein):
		print 'skipping making measuring ISO lOIII5008 as infile does not exist'
		return None
	else:
		# reading file
		hdu=fits.open(filein)
		img=hdu[0].data*u.Unit(hdu[0].header['BUNIT'])
		# setting isophotocut
		isophotocut=isophotocut_base*(1.+obj.sdss.z)**-4

		tabout, imgmos = measureISO(img,isophotocut,smoothing)

		if towriteimgmos:
			fileout=obj.dir_obj+'stamp-lOIII5008_I_mos'+str(smoothing)+'.fits'
			header=hdu[0].header
			headerout=header.copy()
			for col in ['CRPIX1','CRPIX2',]: headerout[col]=header[col]/smoothing
			for col in ['CD1_1','CD1_2','CD2_1','CD2_2',]: headerout[col]=header[col]*smoothing
			prihdu = fits.PrimaryHDU(imgmos, header=headerout)
			prihdu.writeto(fileout, clobber=True)
		return tabout




def measureISO(img,isophotocut,smoothing=4,pixelsize=0.396):
	"""
	PURPOSE: Given isophotocut, measure the nebular size, area, and flux. 
	         To increase singal to noise ratio, smoothing can be applied (recommended). 


	PARAMETERS: 
			img (2d array)     :  line intensity map [erg/s/cm2/arcsec2]
			isophotocut	(float):  isophoto cut [erg/s/cm2/arcsec2]
			smoothing=4:       :  the size of the smoothing superpixels [pix]
			pixelsize=0.396    :  the image pixel size [arcsec]

	DESCRIPTION: 
		 If smoothing > 1, mosaicing (average) each block of smoothing * smoothing pixels 
		 in the original image to one superpixel in the mosaiced image. All measurements are done in
		 the mosaiced image, which has worse resolution but better SNR. 

		 dISO is the largest distance between any pair of superpixel that has flux > isophotocut.
		 rISO is the largest distance between image center and superpixels that has flux > isophotocut.
		 area is the sum of area in the superpixels that has flux > isophotocut.
		 flux is the sum of flux in the superpixels that has flux > isophotocut. (This may not capture all the faint emission, but will not suffer from negative over subtractions)
		
	OUTPUT: 
		tabout, imgmos 

			tabout: a table with columns: dISO, rISO, area, flux
			imgmos: mosaiced and isophoto clipped image
	"""
	# sanity check of unit
	try: img.unit
	except: 
		print "assuming image unit of [erg s-1 cm-2 arcsec-2]"
		imgunit=u.Unit('erg s-1 cm-2 arcsec-2')
	else:  imgunit=img.unit

	# setting
	superpixelsize=pixelsize*smoothing*u.arcsec

	#==== smoothing
	imgmos=mosaicImg(img,smoothing)
	nx,ny=imgmos.shape
	# set all the faint featuers 0, so only bright super pixels will have value
	imgmos[imgmos*imgunit<isophotocut]=0.
	nxm,nym=imgmos.shape
	nbp=sum(imgmos>0.) # count number of bright super pixels

	#==== calculation
	if nbp ==0:
		dISO=0.*superpixelsize
		rISO=0.*superpixelsize
		area=0.*superpixelsize**2
		flux=0.*imgunit*superpixelsize**2
	elif nbp >0:
		# calculate dISO
		ixs,iys=np.where(imgmos>0)
		distances=np.zeros(nbp*nbp)
		for i in range(nbp):
			for j in range(nbp):
				distances[i*nbp+j]=sqrt((ixs[i]-ixs[j])**2+(iys[i]-iys[j])**2)
		dISO=max(distances)*superpixelsize
		# calculate riso
		distances=np.zeros(nbp)
		for i in range(nbp):
			distances[i]=sqrt((ixs[i]-0.5*(nxm-1))**2+(iys[i]-0.5*(nym-1))**2)
		rISO=max(distances)*superpixelsize
		# calculate area
		area=nbp*superpixelsize**2
		# calculate flux
		flux=sum(imgmos)*imgunit*superpixelsize**2

	# output
	tabout=Table()
	for var in ['dISO','rISO','area','flux']: 
		tabout[var]=[eval(var).value]
		tabout[var].unit=eval(var).unit

	return tabout, imgmos


def mosaicImg(img,smoothing):
	"""
	mosaicing image 

	DESCRIPTION: 
	 If smoothing > 1, mosaicing (average) each block of smoothing * smoothing pixels 
	 in the original image to one superpixel in the mosaiced image. 

	"""
	# smoothing
	nx,ny=img.shape
	nxm,nym=int(nx/smoothing),int(ny/smoothing)
	imgmos=np.zeros([nxm,nym])
	for i in range(nxm):
		for j in range(nym):
			imgmos[i,j]=np.average(img[smoothing*i:smoothing*i+smoothing,smoothing*j:smoothing*j+smoothing])
	return imgmos
