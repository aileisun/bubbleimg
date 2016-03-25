# main_mullaney.py
# ALS 2015/08/14

"""
 For selected mullaney targets automatically make the stamp images, color images, and OIII map

"""
import numpy as np
from astropy.table import Table, join

import os

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/')
import catalogue_util

import run_batch
reload(run_batch)


def main():

	#==== read input
	# Sample Mullaney_T2vLOIII.fits
	# containing 360 very Luminous (LOIII>1.e42) Type 2 of all z
	filelistin='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/Mullaney_T2vLOIII.fits'
	listin=Table.read(filelistin,format='fits')

	# #==== set up g - i
	# bandline='g'
	# bandconti='i'

	# blue_edge='OIII'
	# red_edge='Ha'

	# run_mullaney_batch(listin=listin,bandline=bandline,bandconti=bandconti,blue_edge=blue_edge,red_edge=red_edge,list_exclude=[])


	# #==== set up g - r
	# bandline='g'
	# bandconti='r'

	# blue_edge='Ha'
	# red_edge='OIII'

	# run_mullaney_batch(listin=listin,bandline=bandline,bandconti=bandconti,blue_edge=blue_edge,red_edge=red_edge,list_exclude=[])


	#==== set up r - z
	# running band r - z
	# redshift range 0.134609937514 - 0.24674420522
	# Object RA 231.12016 DEC 40.133236 excluded
	# running 194 objects

	bandline='r'
	bandconti='z'

	blue_edge='OIII'
	red_edge='Ha'

	# list_exclude=[796]

	list_exclude=Table()
	list_exclude['RA']=[231.12016]
	list_exclude['DEC']=[40.133236]

	run_mullaney_batch(listin=listin,bandline=bandline,bandconti=bandconti,blue_edge=blue_edge,red_edge=red_edge,list_exclude=list_exclude)


	#==== set up r - i
	# z1 # 0.25036634792500001
	# z2 # 0.34338609970799999

	bandline='r'
	bandconti='i'

	blue_edge='Ha'
	red_edge='OIII'


	run_mullaney_batch(listin=listin,bandline=bandline,bandconti=bandconti,blue_edge=blue_edge,red_edge=red_edge,list_exclude=[])




def get_zrange(bandline='r',bandconti='i',blue_edge='Ha',red_edge='OIII'):
	"""
	return the redshift range given the band set ups

	blue_edge: which line determines the low z cut
	red_edge: which line determines the high z cut
	"""
	zrangesOIII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/OIIIredshiftrange0.6.txt',format='ascii')
	zrangesHaNII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/HaNIIredshiftrange0.2.txt',format='ascii')

	# set low z cut
	if blue_edge=='Ha':
		z1=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z2']
	elif blue_edge=='OIII':
		z1=zrangesOIII[zrangesOIII['band']==bandline][0]['z1']
	else: 
		raise NameError('edge instruction not understood')

	# set high z cut
	if red_edge=='Ha':
		z2=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z1']
	elif red_edge=='OIII':
		z2=zrangesOIII[zrangesOIII['band']==bandline][0]['z2']
	else: 
		raise NameError('edge instruction not understood')

	z1=max(0,z1)
	z2=max(0,z2)

	return z1, z2


def run_mullaney_batch(listin,bandline='r',bandconti='i',blue_edge='Ha',red_edge='OIII',list_exclude=[]):
	"""
	run the mullaney batch of the z range such that 
		bandline contains OIII lines
		and bandconti contains no lines

		blue_edge: which line determines the low z cut ('Ha' 'OIII')
		red_edge: which line determines the high z cut

	"""
	print 'running band '+bandline+' - '+bandconti
	#===== sanity check
	bands=['u','g','r','i','z']
	if (bandline not in bands):
		raise NameError('bandline not recognized')
	if (bandconti not in bands):
		raise NameError('bandconti not recognized')
	if bandline == bandconti:
		raise NameError('bandline and bandconti coincide')

	#===== setting
	batch='mullaney_'+bandline+bandconti+'_T2'
	dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
	if not os.path.isdir(dir_data): os.mkdir(dir_data)

	#==== selection

	z1, z2 = get_zrange(bandline=bandline,bandconti=bandconti,blue_edge=blue_edge,red_edge=red_edge)
	print 'redshift range '+str(z1)+' - '+str(z2)

	selection=(listin['Z']>z1)&(listin['Z']<z2)&(listin['AGN_TYPE']==2)

	if len(list_exclude)>0:
		id_exclude=np.where(catalogue_util.selectionRADEC(listin,list_exclude))[0]
		for i in id_exclude:
			selection[i]=False 
			print 'Object RA '+str(listin['RA'][i])+' DEC '+str(listin['DEC'][i])+' excluded'

	# determine to run or not
	nobj=sum(selection)
	if nobj>0:
		#==== run
		print 'running '+str(nobj)+' objects'
		list_torun=listin[selection]
		run_batch.run_batch(list_torun,dir_data,batch=batch,bandline=bandline,bandconti=bandconti)
	else:
		#==== run
		print 'skipping as no object was found'




# def run_mullaney_ri_T2():
# 	"""

# 	use r - i to bring out the [OIII] lines 

# 	this batch contains 64 objects where
# 	 	OIII falls in r (filter reponse fucntion = 0.6 )  
# 		HaNII is in z (filter reponse fucntion = 0.2). 

# 	in z range
# 				0.25036  		< z <		 0.34338 
# 	(HaNII at the red side of i)	(OIII at the red side of r)

# 	!!! WARNING: Mullaney_allvLOIIIr only go up to z=0.34319901 !!!
# 	"""

# 	#==== set up
# 	bandline='r'
# 	bandconti='i'
# 	batch='mullaney_ri_T2'
# 	# set output
# 	dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
# 	if not os.path.isdir(dir_data): os.mkdir(dir_data)

# 	#==== read input
# 	filelistin='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/Mullaney_allvLOIIIr.fits'
# 	listm=Table.read(filelistin,format='fits')

# 	#==== selection
# 	zrangesOIII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/OIIIredshiftrange0.6.txt',format='ascii')
# 	zrangesHaNII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/HaNIIredshiftrange0.2.txt',format='ascii')

# 	z1=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z2'] # 0.25036634792500001
# 	z2=zrangesOIII[zrangesOIII['band']==bandline][0]['z2'] # 0.34338609970799999

# 	selection=[(listm['Z']>z1)&(listm['Z']<z2)&(listm['AGN_TYPE']==2)] # such that Halpha does not enter z band and T2
# 	sum(selection) # 64 passes the selection
# 	list_torun=listm[selection]

# 	#==== run
# 	run_batch.run_batch(list_torun,dir_data,batch=batch,bandline=bandline,bandconti=bandconti)


# def run_mullaney_rz_T2():
# 	"""
# 	use r - z to bring out the [OIII] lines 

# 	this batch contains 194 objects where
# 	 	OIII falls in r (filter reponse fucntion = 0.6 )  
# 		HaNII is in i (filter reponse fucntion = 0.2). 

# 	in z range
# 				0.1346  		< z <		 0.2467 
# 	(HaNII at the red side of i)	(OIII at the red side of r)

# 	"""

# 	#==== set up
# 	bandline='r'
# 	bandconti='z'
# 	# set output
# 	dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
# 	if not os.path.isdir(dir_data): os.mkdir(dir_data)
# 	batch='mullaney_rz_T2'

# 	#==== read input
# 	filelistin='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/Mullaney_allvLOIIIr.fits'
# 	listm=Table.read(filelistin,format='fits')

# 	#==== selection
# 	zrangesOIII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/OIIIredshiftrange0.6.txt',format='ascii')
# 	zrangesHaNII=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/HaNIIredshiftrange0.2.txt',format='ascii')

# 	z1=zrangesOIII[zrangesOIII['band']==bandline][0]['z1'] # 0.13460993751399999
# 	z2=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z1'] # 0.24674420522000001
# 	selection=(listm['Z']>z1)&(listm['Z']<z2)&(listm['AGN_TYPE']==2) 

# 	# exception
# 	print "WARNING: Exception, object SDSSJ1524+4007 excluded from run due to error"
# 	selection[796]=False 
# 	# This object SDSSJ1524+4007 RA 231.12016, DEC 40.133236
# 	# spSpec-52765-1293-417.fit
# 	# is in Mullaney but not SDSS Image!!! Why is that???
# 	sum(selection) # 194 passes the selection
# 	list_torun=listm[selection]

# 	#==== run
# 	run_batch.run_batch(list_torun,dir_data,batch=batch,bandline=bandline,bandconti=bandconti)

# # run_batch.write_allImg(list_torun,catalog,batch,bandline,bandconti)

# 	# run_batch.write_mags(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix)
# 	# run_batch.write_measureISO_lOIII5008(filelist=filelist,dir_data=dir_batch,catalog=catalog,batch=batch,suffix=suffix)
# 	# # join the measure nebula table with mullaney big table
# 	# if tojoinmullaney: run_batch.joinmullaney(dir_batch)


if __name__ == '__main__':
	main()

