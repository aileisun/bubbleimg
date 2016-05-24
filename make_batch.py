# make_batch.py
# ALS 2015/08/31

"""
Purpose: do all image operations and measurements on a batch of SDSS objects. 

"""


import os
import numpy as np
from astropy.table import Table, vstack #, join, hstack, vstack
import astropy.units as u

# import sys
# sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/')
# import catalogue_util


import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import blobmaps
reload(blobmaps)


# import measurenebula
# reload(measurenebula)

# reload(fromspec)



def make_batch(dir_batch, list_torun, bandline='r', bandconti='z', catalog='mullaney'):
	""" 
	Establish batch directory that contains summary, list, and object 
	directories with line maps.

	operations in order:
	0. make a batch directory dir_data+batch+'/'
	1. write a summary of the batch
	2. write a list of all objects
	3. run batch_makeblobmaps, which load images, spec, and make [OIII] map

	Parameters
	------
	dir_batch: string
		e.g. 'test/', the directory of the batch
	list_torun: list)
		a list of RA, DEC for objects to run through 
	bandline='r'
	bandconti='z'
	catalog='mullaney'

	Write Output
	------
	bunch of files and directories under dir_batch=dir_data+batch+'/', 
	"""
	# set up output directory
	# dir_batch=dir_data+batch+'/'
	if not os.path.isdir(dir_batch): os.mkdir(dir_batch)

	# write list and summary
	batch_writelist(dir_batch, list_torun) 
	batch_writesummary(dir_batch, list_torun, bandline, bandconti, catalog)

	# load IDs
	# batch_loadIDs(dir_batch, list_torun, catalog)

	# load images and map line maps
	batch_makeblobmaps(dir_batch, list_torun, catalog, bandline, bandconti, update=False)


def batch_writelist(dir_batch, list_torun):
	""" Just write down the list of the objects we ran through """
	# set output
	filename=dir_batch+'list.txt'
	# make content
	list_towrite=list_torun['RA','DEC','SDSSNAME']
	list_towrite.rename_column('SDSSNAME','OBJNAME')
	list_towrite=list_towrite['OBJNAME','RA','DEC']
	# write
	list_towrite.write(filename,format='ascii.fixed_width',delimiter='')
	return list_towrite


def batch_writesummary(dir_batch, list_torun, catalog, bandline, bandconti):
	""" write summary of the batch """
	# set output
	filename=dir_batch+'summary.txt'
	# make content
	batch=os.path.basename(os.path.normpath(dir_batch))
	nobj=len(list_torun)
	z_min=list_torun['Z'].min()
	z_max=list_torun['Z'].max()
	# L_15_min=list_torun['WISE_nuLnu_rf15um'].min()
	# L_15_max=list_torun['WISE_nuLnu_rf15um'].max()
	tsum = Table([[catalog],[batch],[bandline],[bandconti],[nobj],[z_min],[z_max]],
	names=['catalog','batch','band_line','bandconti','n_obj','z_min','z_max'])
	# tsum = Table([[catalog],[batch],[bandline],[bandconti],[nobj],[z_min],[z_max],[L_15_min],[L_15_max]],
	# 	names=['catalog','batch','band_line','bandconti','n_obj','z_min','z_max','L_15_min','L_15_max'])
	# write
	tsum.write(filename,format='ascii.fixed_width',delimiter='')


def batch_loadIDs(dir_batch, list_torun, catalog):
	""" 
	load default files for each obj in the list, by calling obsobj. 
	Files: xid.csv, PhotoObj.csv
	"""
	for i in range(len(list_torun)):
		obj=obsobj(listin=list_torun[i],catalog=catalog,dir_parent=dir_batch,towriteID=True)

		print obj.sdssname


def batch_makeblobmaps(dir_batch, list_torun, catalog, bandline,  bandconti, update=False):
	"""
	PURPOSE: run all blobmaps steps to make all images
	"""
	for i in range(len(list_torun)):
		obj=obsobj(listin=list_torun[i], catalog=catalog, dir_parent=dir_batch, towriteID=True)

		if obj.sdss.xid is not None:
			print 'making blob map '+obj.sdssname
			blobmaps.obj_makeblobmaps(obj, bandline=bandline,bandconti=bandconti,update=update)
		else:
			print 'skipping '+obj.sdssname
			os.rmdir(dir_batch+obj.sdssname)
			append_list_exclude(dir_batch, obj.sdssname, obj.ra, obj.dec)


def	append_list_exclude(dir_batch, sdssname, ra, dec):
	"""add an object to the list of exclude"""
	filename=dir_batch+'list_exclude.txt'
 
 	batch=os.path.basename(os.path.normpath(dir_batch))

	tabnew=Table([[batch],[sdssname],[ra],[dec]],names=['batch','OBJNAME','RA','DEC'])

	if os.path.isfile(filename):
		tab=Table.read(filename,format='ascii')
		if tabnew[0] not in tab:
			tabout=vstack([tab,tabnew])
			tabout.write(filename,format='ascii.fixed_width',delimiter='')
	else:
		tabnew.write(filename,format='ascii.fixed_width',delimiter='')





# measurements - continuum magnitudes
# write_batch_contimags(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix)


# # measurements - nebula properties
# write_measureISO_lOIII5008(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix,isophotocut_base=3e-15*u.Unit('erg s-1 cm-2 arcsec-2'),smoothing=4)
# write_measureISO_lOIII5008(filelist,dir_data,catalog,batch,suffix=suffix,isophotocut_base=5.3e-15*u.Unit('erg s-1 cm-2 arcsec-2'),smoothing=3)

# # join the measure nebula table with mullaney big table
# if tojoinmullaney: 
# 	# joinmullaney(dir_batch)
# 	joinmullaney(dir_batch,filenameout='join_ISO_I3e-15_b4',filenamein_mags='mags.fits',filenamein_ISO='measureISO_I3e-15_b4.fits')
# 	joinmullaney(dir_batch,filenameout='join_ISO_I5.3e-15_b3',filenamein_mags='mags.fits',filenamein_ISO='measureISO_I5.3e-15_b3.fits')


# def write_batch_contimags(filelist,dir_data,catalog,batch,suffix=''):
# 	"""
# 	PURPOSE: write table to store the continuum AB magnitudes of objects in list 
# 	"""
# 	from blobmaps import fromspec

# 	fileout=dir_data+'mags'+suffix

# 	if not os.path.isfile(fileout+'.fits'):
# 		# input
# 		listin=Table.read(filelist,format='ascii')
# 		# operation
# 		tabout=Table()
# 		for i in range(len(listin)):
# 			# declare obj
# 			obj=class_obsobj.obsobj(listin[i],catalog=catalog,batch=batch)
# 			print obj.sdssname

# 			# make table
# 			cols=['fiberMag_u','fiberMag_g','fiberMag_r','fiberMag_i','fiberMag_z']
# 			tabfiber=obj.sdss.photoobj[cols][0]
# 			tabsconti=fromspec.obj_measure_contiABmags(obj)
# 			newrow=hstack([Table(listin[i]),tabfiber,tabsconti])
# 			tabout=vstack([tabout,newrow])

# 		# writeout
# 		tabout.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')
# 		tabout.write(fileout+'.fits',format='fits',overwrite=True)

	

# def write_measureISO_lOIII5008(filelist,dir_data,catalog,batch,suffix='',isophotocut_base=3.e-15*u.Unit('erg s-1 cm-2 arcsec-2'),
# 	smoothing=4):
# 	"""
# 	PURPOSE: make a table to measure ISO properties of lOIII5008 line for listed objects

# 	PARAMETERS: 
# 			dir_data (string)  : path to data directory
# 			catalog='mullaney' : 
# 			filelist='list_Mullaney_z_0.135_0.236_T2.txt'
# 			suffix=''
# 			isophotocut_base=3.e-15*u.Unit('erg s-1 cm-2 arcsec-2')
# 			smoothing=4

# 	(could be improved by using list_torun instead of filelist as input)
# 	"""	
# 	# output
# 	fileout=dir_data+'measureISO_I'+str(isophotocut_base.value)+'_b'+str(smoothing)+suffix
# 	# input
# 	listin=Table.read(filelist,format='ascii')

# 	# operation
# 	tabISO=Table()
# 	for i in range(len(listin)):
# 		# declare obj
# 		obj=class_obsobj.obsobj(listin[i],catalog=catalog,batch=batch)
# 		print obj.sdssname
# 		newrow=hstack([Table(listin[i]),Table(measurenebula.obj_measureISO_lOIII5008(obj,isophotocut_base=isophotocut_base,smoothing=smoothing))])
# 		tabISO=vstack([tabISO,newrow])

# 	# writeout
# 	tabISO.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')
# 	tabISO.write(fileout+'.fits',format='fits',overwrite=True)


# def joinmullaney(dir_data,filenameout=None,filenamein_mags='mags.fits',filenamein_ISO='measureISO_I3e-15_b4.fits',filemullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/Mullaney_allvLOIIIr.fits'):
# 	"""
# 	join the table with mullaney table
# 	"""

# 	# setting - get filenameout without extension
# 	if filenameout==None: filenameout='join'

# 	# files to join

# 	tabin_mags=Table.read(dir_data+filenamein_mags,format='fits')
# 	tabin_ISO=Table.read(dir_data+filenamein_ISO,format='fits')
# 	tabin=join(tabin_mags, tabin_ISO,keys=['RA','DEC'])

# 	# operation
# 	tabmullaney=Table.read(filemullaney,format='fits')
# 	del tabmullaney['NAME']
# 	tout=join(tabin,tabmullaney,keys=['RA','DEC'], join_type='left')

# 	# write out
# 	tout.write(dir_data+filenameout+'.fits',format='fits',overwrite=True)
# 	tout.write(dir_data+filenameout+'.csv',format='ascii.csv')

