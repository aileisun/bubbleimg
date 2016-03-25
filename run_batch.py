# run_batch.py
# ALS 2015/08/31

"""
Purpose: do all image operations and measurements on a batch of SDSS objects. 

"""


import os
import numpy as np
from astropy.table import Table, join, hstack, vstack
import astropy.units as u

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/')
import catalogue_util


import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import alignstamp
reload(alignstamp)

import subtractimg
reload(subtractimg)

import imagedisp_util
reload(imagedisp_util)

import measurenebula
reload(measurenebula)

import fromspec
reload(fromspec)



def run_batch(list_torun,dir_data, bandline='r',bandconti='z',batch='',catalog='mullaney',tojoinmullaney=True,suffix=''):
	"""
	PURPOSE: automatically run through a list of objects (batch) to do the following things:

			 1. down load SDSS 5 band frame images, make stamp images
			 2. and make gri color HumVI image. 
			 3. subtract images to make [OIII] map
			 4. measure [OIII] map and compile file dir_batch+"measureISO_3e-15_lOIII5008.txt"
			 5. (optional) join the [OIII] measurement file with the mullaney big table

	PARAMS: 
			list_torun (list): a list of RA, DEC for objects to run through 
			dir_data (string): e.g., dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
			batch='' (string): e.g., 'mullaney_z_0.135_0.236_T2' 
			catalog='mullaney'
			tojoinmullaney=True

	WRITE OUTPUT: 
			bunch of files and directories under dir_batch=dir_data+batch+'/', 
				e.g., '/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/mullaney_z_0.135_0.236_T2/'
	"""
	# set up output directory
	dir_batch=dir_data+batch+'/'
	if not os.path.isdir(dir_batch): os.mkdir(dir_batch)

	# write list and summary
	filelist=dir_batch+'list.txt'
	writelist(list_torun,filelist)
	filesummary=dir_batch+'summary.txt'
	writesummary(filesummary,list_torun,bandline,bandconti,batch,catalog)

	write_allImg(list_torun,catalog,batch,bandline,bandconti)
	# measurements - continuum magnitudes
	write_mags(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix)
	# measurements - nebula properties
	write_measureISO_lOIII5008(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix,isophotocut_base=3e-15*u.Unit('erg s-1 cm-2 arcsec-2',smoothing=4)
	write_measureISO_lOIII5008(filelist,dir_data,catalog,batch,suffix=suffix,isophotocut_base=5.3e-15*u.Unit('erg s-1 cm-2 arcsec-2',smoothing=3)

	# join the measure nebula table with mullaney big table
	if tojoinmullaney: joinmullaney(dir_batch)

def writesummary(filenameout,list_torun,bandline,bandconti,batch,catalog):
	""" write summary of the batch """

	nobj=len(list_torun)
	z_min=list_torun['Z'].min()
	z_max=list_torun['Z'].max()
	L_15_min=list_torun['WISE_nuLnu_rf15um'].min()
	L_15_max=list_torun['WISE_nuLnu_rf15um'].max()

	tsum = Table([[catalog],[batch],[bandline],[bandconti],[nobj],[z_min],[z_max],[L_15_min],[L_15_max]],
		names=['catalog','batch','band_line','bandconti','n_obj','z_min','z_max','L_15_min','L_15_max'])
	tsum.write(filenameout,format='ascii.fixed_width')


def writelist(list_torun,filenameout):
	""" Just write down the list of the objects we ran through """
	list_towrite=list_torun['RA','DEC','SDSSNAME']
	# catalogue_util.addcolSDSSName(list_towrite)
	list_towrite.rename_column('SDSSNAME','NAME')
	list_towrite=list_towrite['NAME','RA','DEC']
	list_towrite.write(filenameout,format='ascii')
	return list_towrite


def write_allImg(list_torun,catalog,batch,bandline,bandconti):
	"""
	PURPOSE: write all images
	"""
	for i in range(len(list_torun)):
		obj=obsobj(list_torun[i],catalog=catalog,batch=batch)

		print obj.sdssname
		objw_makeallimg(obj, bandline=bandline,bandconti=bandconti)

def objw_makeallimg(obj, bandline='r',bandconti='z'):
	"""
	PURPOSE: for a given obsobj object, make all the images. 
	PARAMS: obj
	WIRTE OUTPUT: obj.dir_obj/...
	"""

	# make stamps		
	alignstamp.objw_makeall(obj)
	# make color images
	imagedisp_util.objw_HumVIgriimages(obj)
	# subtract stamps and make OIII maps
	subtractimg.objw_makeall(obj, bandline=bandline,bandconti=bandconti)


def write_mags(filelist,dir_data,catalog,batch,suffix=''):
	"""
	PURPOSE: write table to store the continuum AB magnitudes of objects in list 
	"""
	fileout=dir_data+'mags'+suffix
	# input
	listin=Table.read(filelist,format='ascii')

	# operation
	tabout=Table()
	for i in range(len(listin)):
		# declare obj
		obj=class_obsobj.obsobj(listin[i],catalog=catalog,batch=batch)
		print obj.sdssname

		# make table
		cols=['fiberMag_u','fiberMag_g','fiberMag_r','fiberMag_i','fiberMag_z']
		tabfiber=obj.sdss.photoobj[cols][0]
		tabsconti=fromspec.obj_measure_contiABmags(obj)
		newrow=hstack([Table(listin[i]),tabfiber,tabsconti])
		tabout=vstack([tabout,newrow])

	# writeout
	tabout.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')
	tabout.write(fileout+'.fits',format='fits',overwrite=True)

	

def write_measureISO_lOIII5008(filelist,dir_data,catalog,batch,suffix='',isophotocut_base=3.e-15*u.Unit('erg s-1 cm-2 arcsec-2'),
	smoothing=4):
	"""
	PURPOSE: make a table to measure ISO properties of lOIII5008 line for listed objects

	PARAMETERS: 
			dir_data (string)  : path to data directory
			catalog='mullaney' : 
			filelist='list_Mullaney_z_0.135_0.236_T2.txt'
			suffix=''
			isophotocut_base=3.e-15*u.Unit('erg s-1 cm-2 arcsec-2')
			smoothing=4

	(could be improved by using list_torun instead of filelist as input)
	"""	
	# output
	fileout=dir_data+'measureISO_I'+str(isophotocut_base.value)+'_b'+str(smoothing)+suffix
	# input
	listin=Table.read(filelist,format='ascii')

	# operation
	tabISO=Table()
	for i in range(len(listin)):
		# declare obj
		obj=class_obsobj.obsobj(listin[i],catalog=catalog,batch=batch)
		print obj.sdssname
		newrow=hstack([Table(listin[i]),Table(measurenebula.obj_measureISO_lOIII5008(obj,isophotocut_base=isophotocut_base,smoothing=smoothing))])
		tabISO=vstack([tabISO,newrow])

	# writeout
	tabISO.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')
	tabISO.write(fileout+'.fits',format='fits',overwrite=True)


def joinmullaney(dir_data,filenameout=None,filemullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/Mullaney_allvLOIIIr.fits'):
	"""
	join the table with mullaney table
	"""

	# setting - get filenameout without extension
	if filenameout==None: filenameout='join_mullaney'

	# files to join
	filenamein1='mags.fits'
	filenamein2='measureISO_3e-15_lOIII5008.fits'

	tabin1=Table.read(dir_data+filenamein1,format='fits')
	tabin2=Table.read(dir_data+filenamein2,format='fits')
	tabin=join(tabin1, tabin2,keys=['RA','DEC'])

	# operation
	tabmullaney=Table.read(filemullaney,format='fits')
	del tabmullaney['NAME']
	tout=join(tabin,tabmullaney,keys=['RA','DEC'], join_type='left')

	# write out
	tout.write(dir_data+filenameout+'.fits',format='fits',overwrite=True)
	tout.write(dir_data+filenameout+'.csv',format='ascii.csv')


# def test():

# 	for i in range(len(list_torun)):
# 		obj=obsobj(list_torun[i],catalog='mullaney',tomatchsdss=False)
# 		print obj.sdssname
# 		subtractimg.objw_lOIII5008_I_png(obj)



