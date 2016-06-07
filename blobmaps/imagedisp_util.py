# imagedisp_util.py
# ALS 2015/08/29
"""
A collection of functions for image displays
"""
import os
from astropy.io import fits
import matplotlib.pylab as plt
import scipy


def objw_HumVIgriimages(obj,update=False):
	"""
	PURPOSE: if file does not already exist, use package HumVI to make gri color images from the stamps for the Magellan targets
	PARAMS: 
			obj
			filename='HumVI_gri.png'

	READ INPUT: 	stamp i,r,g images, e.g., obj.dir_obj+'stamp-i.fits' 
	WRITE OUTPUT: 	obj.dir_obj+filename
	"""
	filename='HumVI_gri.png'

	pathHumVI='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/HumVI/HumVI-master/compose.py'

	dir_obj=obj.dir_obj

	if (not os.path.isfile(dir_obj+filename)) or update:
		commandHumVI=pathHumVI+' -s 1.0,1.1,1.0  -p 1.6,1.6  -o '+dir_obj+filename+' '+dir_obj+'stamp-i.fits '+dir_obj+'stamp-r.fits '+dir_obj+'stamp-g.fits'

		os.system(commandHumVI)



def fits_to_image(filenamein,formats='png',saturatebright=0.5,toshow=False):
	"""
	PURPOSE: make fits image into png image
	PARAMS: 
			 filenamein
			 formats='png' (string)
			 saturatebright= 0.5 (float): a number between 0 and 1, marking relative scale, above which to saturate. 
			 toshow=False (bool)
	READ INPUT: 
			 filenamein
	WRITE OUTPUT
			 filenamein+'.png'
	"""
	# setting
	extension='.'+formats
	if filenamein[-5:]=='.fits': filenameout=filenamein[:-5]+extension
	else: filenameout=filenamein+extension

	if not os.path.isfile(filenamein):
		print "skipping "+filenamein+" as in file does not exist"
	else:
		# read in
		image=fits.getdata(filenamein)
		# check
		if len(image.shape)!=2: raise ValueError("Input fits file is not an 2D image")

		if formats=='png':
			scipy.misc.imsave(filnameout, image)
		else:
			plt.close('all')
			plt.figure()
			f, ax = plt.subplots(1)

			ax.imshow(image,origin='lower',cmap='gist_gray',vmin=0.,vmax=image.max()*saturatebright,interpolation='nearest')

			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
			ax.get_xaxis().set_ticks([])
			ax.get_yaxis().set_ticks([])
			plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
			plt.draw()
			if toshow:	plt.show(block=False)
			plt.savefig(filenameout,bbox_inches='tight')

