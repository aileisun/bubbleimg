from pylab import *

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/sample/Mullaney/pyselect/')
import utilities

from astroquery.sdss import SDSS
import ds9

from astropy.table import Table
import os

import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import align


def main():
	# purpose make images for magellan targets

	#==== set up
	list=Table.read('list_Magellan.txt',format='ascii')
	dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/data/Magellan/'
	
	#==== operation
	i=0

	for i in range(len(list)):
	
		# declare obj
		obj=obsobj(list[i])
		print obj.mag.objid

		# make data directory
		dir_obj=dir_data+'M'+str(obj.mag.objid)+'/'
		obj.dir_obj=dir_obj
		if not os.path.isdir(dir_obj):
			os.makedirs(dir_obj)

		# loading sdss images
		for band in ('g','r','i'):
			filename=dir_obj+'frame-'+band+'.fits'
			if not os.path.isfile(filename):
				im = SDSS.get_images(matches=obj.sdss.xid, band=band)
				im[0].writeto(filename)

		# load sdss photoobj
		filename=dir_obj+'PhotoObj.csv'
		if not os.path.isfile(filename):
			obj.sdss.load_photoobj()
			obj.sdss.photoobj.write(filename,format='ascii.csv')
		else:
			obj.sdss.photoobj=Table.read(filename,format='ascii.csv')


		# make aligned stamps
		filename=dir_obj+'images_stamp.npy'
		#if not os.path.isfile(filename):
		images_stamp=align.getalignedstampImages(obj,bands=('g','r','i'), band_rf='r',savefits=True)
		np.save(filename,images_stamp)
		#else:
			#images_stamp=np.load(filename)


		# display images
		d=ds9.ds9('ac17039e:61908')
		for i in range(3):
			print "band ", i
			d.set_np2arr(images_stamp[i])
			raw_input('Enter to continue...')



		# make color images

		HumVI='/Users/aisun/Documents/Astro/Thesis/bbselection/algorithm/display/HumVI/HumVI-master/compose.py  -s 1.0,1.3,0.8  -p 0.3,1.0  -o '+dir_obj+'gri_test.png '+dir_obj+'stamp-i.fits '+dir_obj+'stamp-r.fits '+dir_obj+'stamp-g.fits'


		os.system(HumVI)






#	for i in range(len(list)):
#	
#		# declare obj
#		obj=obsobj(list[i])
#		print obj.mag.objid
#		
#		# make data directory
#		dir_obj=dir_data+'M'+str(obj.mag.objid)+'/'
#		HumVI='/Users/aisun/Documents/Astro/Thesis/bbselection/algorithm/display/HumVI/HumVI-master/compose.py  -s 1.0,1.3,0.8  -p 0.3,1.0  -o '+dir_obj+'gri_test.png '+dir_obj+'stamp-i.fits '+dir_obj+'stamp-r.fits '+dir_obj+'stamp-g.fits'
#
#
#		os.system(HumVI)
