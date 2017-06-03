# imagedisp_util.py
# ALS 2015/08/29
"""
A collection of functions for image displays
"""
import os
import numpy as np
from astropy.io import fits

from skimage.io import imsave


def fits_to_png(filenamein, filenameout=None, vmin=None, vmax=None):
	"""
	make fits image into png image

	PARAMS
	----------
	filenamein (string)
	filenameout = None (string)
		if set to None then it will be filenamin+'.png'
	vmin = None (float)
		the value to set to 0 (black)

	vmax = None (float)
		the value to set to 255 (white)


	READ INPUT
	----------
	filenamein

	WRITE OUTPUT
	----------
	filenameout or filenamein +'.png'
	"""

	# setting filenameout
	if filenameout is None:
		extension='.png'
		if filenamein[-5:]=='.fits': 
			filenameout=filenamein[:-5]+extension
		else: 
			filenameout=filenamein+extension

	if not os.path.isfile(filenamein):
		print "skipping "+filenamein+" as in file does not exist"
	else:
		# read in
		image=fits.getdata(filenamein)

		# setting vmin and vmax
		if vmin is None:
			vmin = np.min(image)
		if vmax is None:		
			vmax = np.max(image)

		image_scale = np.clip((image-vmin) / (vmax-vmin), 0., 1.).astype('f2')

		image_scale_flip = image_scale[::-1, :]

		imsave(filenameout, image_scale_flip)



# from .. import external_links


# def objw_HumVIgriimages(obj, bands='gri', update=False):
# 	"""
# 	PURPOSE: if file does not already exist, use package HumVI to make gri color images from the stamps for the Magellan targets
# 	PARAMS: 
# 			obj
# 			filename='HumVI_gri.png'

# 	READ INPUT: 	stamp i,r,g images, e.g., obj.dir_obj+'stamp-i.fits' 
# 	WRITE OUTPUT: 	obj.dir_obj+filename
# 	"""
# 	filename = 'HumVI_'+bands+'.png'

# 	pathHumVI = external_links.file_humvi_compose
# 	print pathHumVI

# 	dir_obj = obj.dir_obj

# 	if (not os.path.isfile(dir_obj+filename)) or update:
# 		commandHumVI = pathHumVI+' -s 1.0,1.1,1.0  -p 1.6,1.6  -o '+dir_obj+filename+' '+dir_obj+'stamp-'+bands[2]+'.fits '+dir_obj+'stamp-'+bands[1]+'.fits '+dir_obj+'stamp-'+bands[0]+'.fits'

# 		os.system(commandHumVI)
