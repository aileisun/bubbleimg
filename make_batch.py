# make_batch.py
# ALS 2015/08/31

"""
Purpose: 
for a batch (list of objects) download relevent SDSS metadata, spectrum, and images, and perform operation on images to get OIII map. 

"""


import os
import numpy as np
import astropy.table as at
import astropy.units as u


import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import blobmaps
reload(blobmaps)



def make_batch(dir_batch, list_torun, bandline='r', bandconti='z', catalog='mullaney', tosummarize=True):
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
	# set up output directory (dir_batch=dir_data+batch+'/')
	if not os.path.isdir(dir_batch): 
		os.mkdir(dir_batch)

	batch_loadIDs(dir_batch, list_torun, catalog)

	# write list and summary
	batch_writelist(dir_batch, list_torun) 

	if tosummarize:
		batch_writesummary(dir_batch, list_torun, bandline, bandconti, catalog)

	# load images and map line maps
	batch_makeblobmaps(dir_batch, list_torun, catalog, bandline, bandconti, update=False)


def batch_loadIDs(dir_batch, list_torun, catalog):
	""" 
	load default files for each obj in the list, by calling obsobj. 
	Files: xid.csv, PhotoObj.csv
	"""
	for i in range(len(list_torun)):
		obj=obsobj(listin=list_torun[i],catalog=catalog,dir_parent=dir_batch,towriteID=True)

		print obj.sdssname

def batch_writelist(dir_batch, list_torun):
	""" Just write down the list of the objects we ran through """
	# set output
	filename=dir_batch+'list.txt'
	# make content
	# list_towrite=list_torun['RA','DEC','SDSSNAME']
	if 'OBJNAME' not in list_torun.colnames:
		list_torun.rename_column('SDSSNAME','OBJNAME')
	list_towrite=list_torun['OBJNAME','RA','DEC']
	# write
	list_towrite.write(filename,format='ascii.fixed_width',delimiter='')


def batch_writesummary(dir_batch, list_torun, catalog, bandline, bandconti):
	""" write summary of the batch """
	# set output
	filename=dir_batch+'summary.txt'
	# make content
	batch=os.path.basename(os.path.normpath(dir_batch))
	nobj=len(list_torun)

	if 'Z' in list_torun.colnames:
		z_min=list_torun['Z'].min()
		z_max=list_torun['Z'].max()
		tsum = at.Table([[catalog],[batch],[bandline],[bandconti],[nobj],[z_min],[z_max]], names=['catalog','batch','band_line','bandconti','n_obj','z_min','z_max'])
	else: 
		tsum = at.Table([[catalog],[batch],[bandline],[bandconti],[nobj]], names=['catalog','batch','band_line','bandconti','n_obj'])

	# write
	tsum.write(filename,format='ascii.fixed_width',delimiter='')


def batch_makeblobmaps(dir_batch, list_torun, catalog, bandline,  bandconti, update=False):
	"""
	PURPOSE: run all blobmaps steps to make all images
	"""

	if os.path.isfile(dir_batch+'list_exclude.txt'):
		listexclude = at.Table.read(dir_batch+'list_exclude.txt',format='ascii')
	else: 
		listexclude = at.Table(data=[[]], names=['OBJNAME'])

	for i in range(len(list_torun)):
		objname = list_torun['OBJNAME'][i]
		if objname not in listexclude['OBJNAME']:
			print "calling obsobj"
			obj = obsobj(listin=list_torun[i], catalog=catalog, dir_parent=dir_batch, towriteID=True)

			if obj.sdss.xid is not None:
				print 'making blob map '+obj.sdssname
				blobmaps.obj_makeblobmaps(obj, bandline=bandline,bandconti=bandconti,update=update)
				print 'making scaleconti '+obj.sdssname
				blobmaps.makemap.objw_makecontiscale(obj, bandline=bandline,bandconti=bandconti,update=update)
			else:
				print 'skipping '+obj.sdssname
				os.rmdir(dir_batch+obj.sdssname)
				append_list_exclude(dir_batch, obj.sdssname, obj.ra, obj.dec)


def	append_list_exclude(dir_batch, sdssname, ra, dec):
	"""add an object to the list of exclude"""
	filename=dir_batch+'list_exclude.txt'
 
 	batch=os.path.basename(os.path.normpath(dir_batch))

	tabnew=at.Table([[batch],[sdssname],[ra],[dec]],names=['batch','OBJNAME','RA','DEC'])

	if os.path.isfile(filename):
		tab=at.Table.read(filename,format='ascii')
		if tabnew[0] not in tab:
			tabout=at.vstack([tab,tabnew])
			tabout.write(filename,format='ascii.fixed_width',delimiter='')
	else:
		tabnew.write(filename,format='ascii.fixed_width',delimiter='')

