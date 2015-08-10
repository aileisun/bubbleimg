from pylab import *
from astroquery.sdss import SDSS
import ds9
from astropy.table import Table
import os
import pickle as pickle

import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import align
reload(align)


dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/Magellan/'
filelist='list_Magellan.txt'

pathHumVI='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/HumVI/HumVI-master/compose.py'

def makeStampsImages(stampbands=['u','g','r','i','z']):
	"""
	download frame images and make stamp images for the Magellan targets
	"""

	#==== set up
	listmagellan=Table.read(filelist,format='ascii')
	tods9 =False
	
	#==== operation
	#i=0

	for i in range(len(listmagellan)):
	
		# declare obj
		obj=obsobj(listmagellan[i])
		print obj.galaxy.OBJID
		

		# make data directory
		dir_obj=dir_data+'M'+str(obj.galaxy.OBJID)+'/'
		obj.dir_obj=dir_obj
		if not os.path.isdir(dir_obj):
			os.makedirs(dir_obj)

		# saving sdss frame images
		for band in stampbands:
			filename=dir_obj+'frame-'+band+'.fits'
			if not os.path.isfile(filename):
				im = SDSS.get_images(matches=obj.sdss.xid, band=band)
				im[0].writeto(filename)

		# saving sdss photoobj
		filename=dir_obj+'PhotoObj.csv'
		obj.sdss.photoobj.write(filename,format='ascii',delimiter=',')

		# saving aligned stamps
		filename=dir_obj+'images_stamp.npy'
		images_stamp=align.getalignedstampImages(obj,bands=stampbands, band_rf='r',savefits=True)
		np.save(filename,images_stamp)


		if tods9:
			# display images
			d=ds9.ds9('ac17039e:61908')
			for i in range(len(stampbands)):
				print "band ", stampbands[i]
				d.set_np2arr(images_stamp[i])
				raw_input('Enter to continue...')


def makeStampsImages_whenSDSSqueryfailed(stampbands=['u','g','r','i','z']):
	"""
	PURPOSE: download frame images and make stamp images for the Magellan targets
		if SDSS.query_region() in class_obsobj failed, and xid cannot be retrieved, then use this funciton
	"""
	class obsobj_fake():
		pass
	class sdss():
		pass

	#==== set up
	listmagellan=Table.read(filelist,format='ascii')
	
	for i in range(len(listmagellan)):
	
		# declare obj
		obj=obsobj_fake()
		obj.ra=listmagellan[i]['RA']
		obj.dec=listmagellan[i]['DEC']
		obj.objid=listmagellan[i]['NAME']
		obj.dir_obj=dir_data+listmagellan[i]['NAME']+'/'
		obj.xid=Table.read(obj.dir_obj+'PhotoObj.csv',format='ascii')
		obj.sdss=sdss()
		obj.sdss.photoobj=Table.read(obj.dir_obj+'PhotoObj.csv',format='ascii')
		obj.sdss.ra=obj.sdss.photoobj['ra'].data[0]
		obj.sdss.dec=obj.sdss.photoobj['dec'].data[0]

		print obj.objid
		

		# saving sdss frame images
		for band in stampbands:
			filename=obj.dir_obj+'frame-'+band+'.fits'
			if not os.path.isfile(filename):
				im = SDSS.get_images(matches=obj.xid, band=band)
				im[0].writeto(filename)

		# saving aligned stamps
		filename=obj.dir_obj+'images_stamp.npy'
		images_stamp=align.getalignedstampImages(obj,bands=stampbands, band_rf='r',savefits=True)
		np.save(filename,images_stamp)





def makeHumVIgriimages(filename='HumVI_gri.png'):
	"""
	make HumVI gri color images from the stamps for the Magellan targets
	"""
	listmagellan=Table.read(filelist,format='ascii')

	for i in range(len(listmagellan)):
		dir_obj=dir_data+listmagellan[i]['NAME']+'/'

		commandHumVI=pathHumVI+' -s 1.0,1.3,0.8  -p 0.3,1.0  -o '+dir_obj+filename+' '+dir_obj+'stamp-i.fits '+dir_obj+'stamp-r.fits '+dir_obj+'stamp-g.fits'

		os.system(commandHumVI)

