"""
WARNING: SDSS.query_region failed
"""


from pylab import *
import os
import sqlcl

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/')
import catalogue_util

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from astropy.io import fits
import astropy.units as u

dir_data_SDSS='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/SDSS/'

dir_data_magellan='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/magellan/'
listname_magellan='/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/sample/Magellan1406_list.txt'
list_magellan=Table.read(listname_magellan,format='ascii')

dir_data_mullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
listname_mullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/ALPAKA_v1.fits'
list_mullaney=Table.read(listname_mullaney,format='fits')

class obsobj(object):
	""" 
	DESCRIPTION: Class of SDSS Object identified by ra, dec 
				 It contains SDSS PhotoObj info and other magnitude infos 
				 (Only objects both in photoobj and specobj are recognized)

	SYNTAX: class_obsobj.obsobj(listin,catalog='magellan')

	PARAMETERS: 
			listin (table row): is a table row containing two columns ['RA','DEC'] in [degree].
			catalog (string): name of the catalog, one of the following:
								  "magellan": output directory will be, for example, /SDSS/data/magellan/M2100. 
								  "mullaney": output directory will be, for example, /SDSS/data/mullaney/JXXXX+XXXX. 
	"""
	
	def __init__(self,listin,catalog='magellan',tomatchsdss=True,batch=''):
		# Initialize SDSSObj with ra dec
		if 'ra' in listin.dtype.names:
			self.ra=listin['ra']
			self.dec=listin['dec']
		if 'RA' in listin.dtype.names:
			self.ra=listin['RA']
			self.dec=listin['DEC']

		#== get SDSS Name
		self.sdssname=catalogue_util.getSDSSName_fromlist(listin)
		self.dir_obj=dir_data_SDSS+self.sdssname+'/'

		#== match wit sdss
		if tomatchsdss:
			self.sdss = obsobj.sdss(self)

		#== match wit magellan 
		if catalog=='magellan':
			print "[obsobj] matching with Magellan 1406 sample"
			row=list_magellan[all([absolute(list_magellan['RA']-self.ra)<0.001,absolute(list_magellan['DEC']-self.dec)<0.001],axis=0)]
			if len(row) == 1: 
				self.galaxy=galaxy.galaxy(row['OBJID'].data[0])
				self.dir_obj=dir_data_magellan+'M'+str(self.galaxy.OBJID)+'/'
				print "matched to "+str(self.galaxy.OBJID)
			else: raise NameError('no match or duplicate to Magellan 1406 sample')

		#== match with mullaney
		elif catalog=='mullaney':
			print "[obsobj] matching with Mullaney+13 sample"
			row=catalogue_util.selectRADEC(list_mullaney,[listin],radius=3.,verbos=True)

			if batch!='':batch=batch+'/'
			if len(row) == 1: 
				self.mullaney=row
				self.dir_obj=dir_data_mullaney+batch+self.sdssname+'/'
				print "matched to "+row['NAME'][0]
			else: raise NameError('no match or duplicate to Mullaney 13 sample')



	
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
			self.outer=outerself

		def print_xid(self):
			print self.xid

		def get_images(self):
			return SDSS.get_images(matches=self.xid)

		def get_spectra(self):
			"""
			PURPOSE: astroquery SDSS spectrum
			RETURN OUTPUT: sp (hdulist)
				----- HDUList ---- 
				HDU 0  : Header info from spPlate
			    HDU 1  : Coadded spectrum from spPlate   <--- THIS IS WHAT I NEED
			    		 Columns: ['flux', 'loglam', 'ivar', 'and_mask', 'or_mask', 'wdisp', 'sky', 'model']
			    HDU 2  : Summary metadata copied from spAll
			    HDU 3  : Line fitting metadata from spZline
			    HDU 4+ : [Optional] Individual spCFrame spectra [B, R for each exposure]

				For more informaiton, see http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
			"""
			return SDSS.get_spectra(matches=self.xid)

		def load_spectra(self, towritefits=True):
			"""
			PURPOSE: If self.dir_obj/spec.fits does not exist locally, download and store 
			         SDSS spectrum of the objec. 
					 If spec.fits exists locally, then load spec from local file. 

			EXAMPLES: 
					spectable=obj.sdss.load_spectra(obj)[1].data
					spec, lcoord= spectable['flux'], 10.**spectable['loglam']

			PARAMETERS:		obj (object of obsobs)
							towritefits (bool)
			RETURN OUTPUT: sp (hdulist)
			WRITE OUTPUT:  obj.dir_obj+'spec.fits'

				----- HDUList ---- 
				HDU 0  : Header info from spPlate
			    HDU 1  : Coadded spectrum from spPlate   <--- ((THIS IS WHAT I NEED))
			    		 Columns: ['flux', 'loglam', 'ivar', 'and_mask', 'or_mask', 'wdisp', 'sky', 'model']
			    HDU 2  : Summary metadata copied from spAll
			    HDU 3  : Line fitting metadata from spZline
			    HDU 4+ : [Optional] Individual spCFrame spectra [B, R for each exposure]

				For more informaiton, see http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
			"""
			filename=self.outer.dir_obj+'spec.fits'
			if os.path.isfile(filename): sp=fits.open(filename)
			else: 
				sp = SDSS.get_spectra(matches=self.xid)
				if len(sp)!=1: raise ValueError("SDSS spec obj not uniquely identified. ")
				else: sp=sp[0]
				if towritefits: sp.writeto(filename)
			return sp

		def get_speclcoord(self, wunit=False):
			"""
			PURPOSE: retrieve SDSS spec of and obj as spec and lcoord arrays.

			EXAMPLES: 					
					spec, lcoord = get_spec_lcoord(obj)

			PARAMETERS:		obj (object of obsobs)
							wunit (bool)  to attach unit or not

			RETURN OUTPUT: spec (nparray), lcoord (nparray)
			"""
			spectable=self.load_spectra(self)[1].data
			spec, lcoord= spectable['flux'], 10.**spectable['loglam']

			if not wunit:	
				return spec, lcoord
			else:	
				u_spec=1.e-17*u.Unit('erg / (Angstrom cm2 s)')
				u_lcoord=u.AA
				return spec*u_spec, lcoord*u_lcoord


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
	# 		row=list_magellan[all([absolute(list_magellan['RA']-outerself.ra)<0.001,absolute(list_magellan['DEC']-outerself.dec)<0.001],axis=0)]
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

