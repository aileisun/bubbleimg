# Magellan_main.py
# ALS reorganized 2015/08/10

"""
 For Magellan targets automatically make the stamp images, color images, and extract color along slits

"""


from pylab import *
from astropy.table import Table, hstack, vstack
import os

import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import alignstamp
reload(alignstamp)

import subtractimg
reload(subtractimg)

import fromspec
reload(fromspec)


import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy

import magellan_getslit
reload(magellan_getslit)

import imagedisp_util
reload(imagedisp_util)

import measurenebula
reload(measurenebula)


dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/magellan/'



def main():
	"""
	automatically make the stamp images, color images, and extract color along slits for Magellan targets

	products saved to self.dir_obj=magellan_dir_data+'M'+str(self.galaxy.OBJID)+'/'
	"""

	#==== set up
	filelist=dir_data+'list_Magellan.txt'
	listmagellan=Table.read(filelist,format='ascii')

	for i in range(len(listmagellan)):
		# declare obj
		obj=obsobj(listmagellan[i],catalog='magellan')
		print obj.galaxy.OBJID

		# make stamps		
		alignstamp.objw_makeall(obj)
		# make color images
		imagedisp_util.objw_HumVIgriimages(obj)
		# subtract stamps and make OIII maps
		subtractimg.objw_makeall(obj)

	# Extract info from stamps along Magellan slits
	galaxy.posIterator(magellan_getslit.poswSlitFluxDensity)

	kwargs={'linetag':'lOIII5008'}
	galaxy.posIterator(magellan_getslit.poswSlitLineIntensity,**kwargs)

	galaxy.posIterator(magellan_getslit.poswSlitAll)

	measurenebula.write_measureISO_lOIII5008(dir_data=dir_data,catalog='magellan',filelist='list_Magellan.txt')

	
	# makemeasureISO_lOIII5008()

# def makemeasureISO_lOIII5008():
# 	"""
# 	make a table to measure ISO properties of lOIII5008 line for all magellan objects
# 	"""
# 	import astropy.units as u
# 	# settings
# 	isophotocut_base=3.e-15*u.Unit('erg s-1 cm-2 arcsec-2')

# 	# output
# 	fileout=dir_data+'measureISO_'+str(isophotocut_base.value)+'_lOIII5008'
# 	# input
# 	filelist='list_Magellan.txt'
# 	listmagellan=Table.read(filelist,format='ascii')

# 	# operation
# 	tabISO=Table()
# 	for i in range(len(listmagellan)):
# 		# declare obj
# 		obj=obsobj(listmagellan[i],catalog='magellan')
# 		print obj.galaxy.OBJID

# 		newrow=hstack([Table(listmagellan[i]),Table(measurenebula.obj_measureISO_lOIII5008(obj,isophotocut_base=isophotocut_base))])
# 		tabISO=vstack([tabISO,newrow])

# 	# writeout
# 	tabISO.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')
# 	tabISO.write(fileout+'.fits',format='fits',overwrite=True)


# def makeStampsImages(obj,stampbands=['u','g','r','i','z']):
# 	"""
# 	download frame images and make stamp images for all the Magellan targets
# 	"""
	
# 	alignstamp.alignstamp(obj, bands=stampbands, band_rf='r', xwidth=57, ywidth=57, tods9=False)



# def makeSubtractedImages(obj):
# 	"""
# 	make all the r minus other band subtracted images. 
# 	"""
# 	fromspec.plotSpecwFilters(obj)
# 	subtractimg.objw_stamp_band_minusconti(obj, band1='r', band2='u', savefits=True)
# 	subtractimg.objw_stamp_band_minusconti(obj, band1='r', band2='g', savefits=True)
# 	subtractimg.objw_stamp_band_minusconti(obj, band1='r', band2='i', savefits=True)
# 	subtractimg.objw_stamp_band_minusconti(obj, band1='r', band2='z', savefits=True)
# 	subtractimg.objw_stamp_lOIII5008_F(obj)
# 	subtractimg.objw_stamp_lOIII5008_I(obj)


# #================ Deprecated ===============
# def makeStampsImages_whenSDSSqueryfailed(listmagellan,stampbands=['u','g','r','i','z']):
# 	"""
# 	PURPOSE: download frame images and make stamp images for the Magellan targets
# 		if SDSS.query_region() in class_obsobj failed, and xid cannot be retrieved, then use this funciton
# 	"""
# 	class obsobj_fake():
# 		pass
# 	class sdss():
# 		pass
	
# 	for i in range(len(listmagellan)):
	
# 		# declare obj
# 		obj=obsobj_fake()
# 		obj.ra=listmagellan[i]['RA']
# 		obj.dec=listmagellan[i]['DEC']
# 		obj.objid=listmagellan[i]['NAME']
# 		obj.dir_obj=dir_data+listmagellan[i]['NAME']+'/'
# 		obj.xid=Table.read(obj.dir_obj+'PhotoObj.csv',format='ascii')
# 		obj.sdss=sdss()
# 		obj.sdss.photoobj=Table.read(obj.dir_obj+'PhotoObj.csv',format='ascii')
# 		obj.sdss.ra=obj.sdss.photoobj['ra'].data[0]
# 		obj.sdss.dec=obj.sdss.photoobj['dec'].data[0]

# 		print obj.objid

# 		alignstamp.alignstamp(obj, bands=stampbands, band_rf='r', xwidth=57, ywidth=57, tods9=False)
# 		# # saving sdss frame images
# 		# for band in stampbands:
# 		# 	filename=obj.dir_obj+'frame-'+band+'.fits'
# 		# 	if not os.path.isfile(filename):
# 		# 		im = SDSS.get_images(matches=obj.xid, band=band)
# 		# 		im[0].writeto(filename)

# 		# # saving aligned stamps
# 		# filename=obj.dir_obj+'images_stamp.npy'
# 		# images_stamp=align.getalignedstampImages(obj,bands=stampbands, band_rf='r',savefits=True)
# 		# np.save(filename,images_stamp)

if __name__ == '__main__':
	main()
