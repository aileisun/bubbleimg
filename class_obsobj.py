from pylab import *

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/sample/Mullaney/pyselect/')
import utilities

import sqlcl

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

class obsobj(object):
	""" Class of SDSS Object identified by ra, dec 
		
		Only objects both in photoobj and specobj are recognized
	"""
	
	def __init__(self,list):
		# Initialize SDSSObj with ra dec
		if 'ra' in list.dtype.names:
			self.ra=list['ra']
			self.dec=list['dec']
		if 'RA' in list.dtype.names:
			self.ra=list['RA']
			self.dec=list['DEC']
		self.sdss = obsobj.sdss(self)
		self.mag = obsobj.mag(self)
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

		def print_xid(self):
			print self.xid

		def get_images(self):
			return SDSS.get_images(matches=self.xid)

		def get_spectra(self):
			return SDSS.get_spectra(matches=self.xid)


		def load_photoobj(self):
			# load photoobj table and store it in obj.sdss.photoobj
			sql="SELECT p.* FROM PhotoObj AS p WHERE p.objid="+str(self.objid)
			#txt=sqlcl.query(sqlc).readlines()
			txt=sqlcl.query(sql).read()

			f=open('PhotoObj.txt','w')
			f.write(txt)
			f.close()
			self.photoobj=Table.read('PhotoObj.txt',format='ascii')

			# image position related quantities:
			#sql="SELECT p.colc, p.colcErr, p.colc_u, p.colc_g, p.colc_r, p.colc_i, p.colc_z, p.colcErr_u, p.colcErr_g, p.colcErr_r, p.colcErr_i, p.colcErr_z, p.rowc, p.rowcErr, p.rowc_u, p.rowc_g, p.rowc_r, p.rowc_i, p.rowc_z, p.rowcErr_u, p.rowcErr_g, p.rowcErr_r, p.rowcErr_i, p.rowcErr_z FROM PhotoObj AS p WHERE p.objid="+str(obj.sdss.objid)


	class mag():

		def __init__(self,outerself):
			print "[obsobj] matching with Magellan"
			filename='/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/sample/Magellan1406.fits'
			table=Table.read(filename,format='fits')
			row=table[all([absolute(table['RA']-outerself.ra)<0.001,absolute(table['DEC']-outerself.dec)<0.001],axis=0)]
			if len(row) == 1:
				self.objid=row['OBJID'][0]
				self.AGN_TYPE=row['AGN_TYPE'][0]
				self.SDSSOIII_5007T_LUM=row['OIII_5007T_LUM'][0]
				self.SDSSOIII_5007AVG_FWHM=row['OIII_5007AVG_FWHM'][0]
				self.W1_Fnu=row['W1_Fnu'][0]
				self.W2_Fnu=row['W2_Fnu'][0]
				self.W3_Fnu=row['W3_Fnu'][0]
				self.W4_Fnu=row['W4_Fnu'][0]
				self.W1_nuLnu_rf=row['W1_nuLnu_rf'][0]
				self.W2_nuLnu_rf=row['W2_nuLnu_rf'][0]
				self.W3_nuLnu_rf=row['W3_nuLnu_rf'][0]
				self.W4_nuLnu_rf=row['W4_nuLnu_rf'][0]
				self.WISE_Sep=row['WISE_Sep'][0]
				self.WISE_nsource=row['WISE_nsource'][0]
				self.WISE_nuLnu_rf8um=row['WISE_nuLnu_rf8um'][0]
				self.WISE_nuLnu_rf20um=row['WISE_nuLnu_rf20um'][0]
				#... more instances can be extracted from sample table
				
			else:
				raise NameError('no match to Magellan 1406 sample')

			# To match Magellan measured line data ... currently not working
			#filename='/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/lines/line_sample.fits'
			#tab=Table.read(filename,format='fits')
			#row=tab[where(tab['OBJID']==str(self.OBJID))]
			#if len(row) == 1:
			#	self.NIGHT=row['NIGHT']
			#	self.FRAME=row['FRAME']
			#	self.SLIT=row['SLIT']
			#	self.OIII5008FWHM=row['OIII5008FWHM']
			#	self.OIII5008FW02M=row['OIII5008FW02M']
			#	self.OIII5008FW01M=row['OIII5008FW01M']
			#	self.OIII5008FW005M=row['OIII5008FW005M']
			#	self.OIII5008RISO=row['OIII5008RISO']
			#	self.OIII5008DISO=row['OIII5008DISO']
			#	self.OIII5008L=row['OIII5008L']
			#else:
			#	raise NameError('no match to Magellan 1406 line')



