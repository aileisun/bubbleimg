"""
WARNING: SDSS.query_region failed
"""


from pylab import *

import os


import sqlcl


import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy


from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

filename='/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/sample/Magellan1406_list.txt'
magellan_list=Table.read(filename,format='ascii')


class obsobj(object):
	""" 
	DESCRIPTION: Class of SDSS Object identified by ra, dec 
				 It contains SDSS PhotoObj info and other magnitude infos 
				 (Only objects both in photoobj and specobj are recognized)

	SYNTAX: class_obsobj.obsobj(list)
			list is a table row containing two columns ['RA','DEC'] in degree.
	"""
	
	def __init__(self,list):
		# Initialize SDSSObj with ra dec
		if 'ra' in list.dtype.names:
			self.ra=list['ra']
			self.dec=list['dec']
		if 'RA' in list.dtype.names:
			self.ra=list['RA']
			self.dec=list['DEC']
		#== match wit sdss
		self.sdss = obsobj.sdss(self)

		#== match wit magellan 
		print "[obsobj] matching with Magellan 1406 sample"
		row=magellan_list[all([absolute(magellan_list['RA']-self.ra)<0.001,absolute(magellan_list['DEC']-self.dec)<0.001],axis=0)]
		if len(row) == 1: self.galaxy=galaxy.galaxy(row['OBJID'].data[0])
		else: raise NameError('no match to Magellan 1406 sample')

		# self.magellan = obsobj.magellan(self)

		#self.crossid()



	
	class sdss():
		def __init__(self,outerself):
			print "[obsobj] matching with SDSS"
			# Retrieving  sdss ids of the closest position-matched spectroscopic object
			
			c = SkyCoord(outerself.ra, outerself.dec, 'icrs', unit='deg')
			result= SDSS.query_region(c,spectro=True) # I hope it pickes up sciencePrimary...
			xid=Table(result[len(result)-1])
			
			self.ra = xid['ra'][0]
			self.dec = xid['dec'][0]
			self.xid	=xid
			self.objid  =xid['objid'][0]
			self.run	=xid['run'][0]
			self.rerun  =xid['rerun'][0]
			self.camcol =xid['camcol'][0]
			self.field  =xid['field'][0]
			self.z      =xid['z'][0]
			self.plate  =xid['plate'][0]
			self.mjd    =xid['mjd'][0]
			self.fiberID=xid['fiberID'][0]
			self.specobjid=xid['specobjid'][0]
			self.run2d=xid['run2d'][0]
			self.instrument=xid['instrument'][0]
			self.load_photoobj()

		def print_xid(self):
			print self.xid

		def get_images(self):
			return SDSS.get_images(matches=self.xid)

		def get_spectra(self):
			return SDSS.get_spectra(matches=self.xid)


		def load_photoobj(self):
			fphotoobj='PhotoObj.txt'
			# load photoobj table and store it in obj.sdss.photoobj
			sql="SELECT p.* FROM PhotoObj AS p WHERE p.objid="+str(self.objid)
			#txt=sqlcl.query(sqlc).readlines()
			txt=sqlcl.query(sql).read()

			f=open(fphotoobj,'w')
			f.write(txt)
			f.close()
			self.photoobj=Table.read(fphotoobj,format='ascii')
			os.remove(fphotoobj)

			# image position related quantities:
			#sql="SELECT p.colc, p.colcErr, p.colc_u, p.colc_g, p.colc_r, p.colc_i, p.colc_z, p.colcErr_u, p.colcErr_g, p.colcErr_r, p.colcErr_i, p.colcErr_z, p.rowc, p.rowcErr, p.rowc_u, p.rowc_g, p.rowc_r, p.rowc_i, p.rowc_z, p.rowcErr_u, p.rowcErr_g, p.rowcErr_r, p.rowcErr_i, p.rowcErr_z FROM PhotoObj AS p WHERE p.objid="+str(obj.sdss.objid)


	# # =========== the class to match with Magellan - currently under construction
	# class magellan():

	# 	def __init__(self,outerself):
	# 		print "[obsobj] matching with Magellan"
	# 		row=magellan_list[all([absolute(magellan_list['RA']-outerself.ra)<0.001,absolute(magellan_list['DEC']-outerself.dec)<0.001],axis=0)]
	# 		if len(row) == 1:
	# 			self.objid=row['OBJID'][0]
	# 			self.galaxy=galaxy.galaxy(self.objid)
	# 			# self.AGN_TYPE=row['AGN_TYPE'][0]
	# 			# self.SDSSOIII_5007T_LUM=row['OIII_5007T_LUM'][0]
	# 			# self.SDSSOIII_5007AVG_FWHM=row['OIII_5007AVG_FWHM'][0]
	# 			# self.W1_Fnu=row['W1_Fnu'][0]
	# 			# self.W2_Fnu=row['W2_Fnu'][0]
	# 			# self.W3_Fnu=row['W3_Fnu'][0]
	# 			# self.W4_Fnu=row['W4_Fnu'][0]
	# 			# self.W1_nuLnu_rf=row['W1_nuLnu_rf'][0]
	# 			# self.W2_nuLnu_rf=row['W2_nuLnu_rf'][0]
	# 			# self.W3_nuLnu_rf=row['W3_nuLnu_rf'][0]
	# 			# self.W4_nuLnu_rf=row['W4_nuLnu_rf'][0]
	# 			# self.WISE_Sep=row['WISE_Sep'][0]
	# 			# self.WISE_nsource=row['WISE_nsource'][0]
	# 			# self.WISE_nuLnu_rf8um=row['WISE_nuLnu_rf8um'][0]
	# 			# self.WISE_nuLnu_rf20um=row['WISE_nuLnu_rf20um'][0]
	# 			# #... more instances can be extracted from sample table
				
	# 		else:
	# 			raise NameError('no match to Magellan 1406 sample')

