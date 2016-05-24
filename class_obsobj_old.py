"""
WARNING: SDSS.query_region failed
"""


from pylab import *
import os

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from astropy.io import fits
import astropy.units as u
import sqlcl

import external_links
import sys
sys.path.append(external_links.pack_magellan_analysis)
import galaxy

import sys
sys.path.append(external_links.pack_catalogue)
import catalogue_util
reload(catalogue_util)
# dir_data_SDSS='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/SDSS/'
# dir_data_magellan='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/magellan/'
# file_list_magellan='/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/sample/Magellan1406_list.txt'
# dir_data_mullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/'
# file_list_mullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/ALPAKA_v1.fits'

list_magellan=Table.read(external_links.file_list_magellan,format='ascii')

list_mullaney=Table.read(external_links.file_list_mullaney,format='fits')

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
		# dir_obj could be overwrite by magellan or mullaney
		self.dir_obj=external_links.dir_data_SDSS+self.sdssname+'/'

		#== match wit magellan 
		if catalog=='magellan':
			print "[obsobj] matching with Magellan 1406 sample"
			row=list_magellan[all([absolute(list_magellan['RA']-self.ra)<0.001,absolute(list_magellan['DEC']-self.dec)<0.001],axis=0)]
			if len(row) == 1: 
				self.galaxy=galaxy.galaxy(row['OBJID'].data[0])
				self.dir_obj=external_links.dir_data_magellan+'M'+str(self.galaxy.OBJID)+'/'
				print "matched to "+str(self.galaxy.OBJID)
			else: raise NameError('no match or duplicate to Magellan 1406 sample')

		#== match with mullaney
		elif catalog=='mullaney':
			print "[obsobj] matching with Mullaney+13 sample"
			row=catalogue_util.selectRADEC(list_mullaney,[listin],radius=3.,verbos=True)

			if batch!='':batch=batch+'/'
			if len(row) == 1: 
				self.mullaney=row
				self.dir_obj=external_links.dir_data_mullaney+batch+self.sdssname+'/'
				print "matched to "+row['NAME'][0]
			else: raise NameError('no match or duplicate to Mullaney 13 sample')

		#== match wit sdss
		if tomatchsdss:
			self.sdss = obsobj.sdss(self)



	
	class sdss():
		def __init__(self,outerself):
			"""
			Initialize object sdss and assign attributes:
				self.outer: outer object, that contains outer.dir_obj
				self.photoboj: table
				self.xid: table 
					columns:
					['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field', 
					'z', 'plate', 'mjd', 'fiberID', 'specobjid', 'run2d', 
					'instrument']
				self.ra: float
				self.dec: float
					..., etc xid parameters.
			"""
			self.outer=outerself
			self.load_xid()
			self.load_photoobj()

		def load_xid(self,tosavexid=True):
			"""
			assign xid (table) an its columns to self as attributes. 
			see get_xid
			"""
			xid=self.get_xid(tosavexid=tosavexid)

			for col in xid.colnames: # setting xid attributes 
				setattr(self, col, xid[col][0])

			self.xid=xid

		def get_xid(self,tosavexid=True):
			"""
			return xid.
			Read xid locally if self.outer.dir_obj+'xid.csv' exist. Otherwise 
			query sdss. 

			Parameters
			------
			self: obj
				contains: 
				self.outer.dir_obj, self.outer.ra, self.outer.dec
			tosavexid=True

			Returns:
			------
			xid: table

			Write Output: (optional)
			------
			self.outer.dir_obj+'xid.csv'
			"""
			filename=self.outer.dir_obj+'xid.csv'

			if os.path.isfile(filename): # retrieve xid locally
				print "[obsobj] reading xid locally"
				xid=Table.read(filename,format='ascii.csv',comment='#')
			else: # download xid from sdss
				print "[obsobj] querying xid from SDSS"
				c = SkyCoord(self.outer.ra, self.outer.dec, 'icrs', unit='deg')
				# I hope it pickes up sciencePrimary...
				result= SDSS.query_region(c,spectro=True) 
				# Retrieving  sdss ids of the closest position-matched spectroscopic object
				xid=Table(result[len(result)-1])
				print "[obsobj] object found"
				if tosavexid and os.path.isdir(os.path.dirname(filename)):
					xid.write(filename,format='ascii.csv')
			return xid


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
			If self.dir_obj/spec.fits exist, load spec locally. Otherwise
			download from SDSS and save it locally. 

			Examples
			------
			spectable=obj.sdss.load_spectra(obj)[1].data
			spec, lcoord= spectable['flux'], 10.**spectable['loglam']

			Parameters
			------
			obj (object of obsobs)
			towritefits (bool)

			Returns
			------
			sp (hdulist)

			Write Outputs
			------
			obj.dir_obj+'spec.fits'
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


		def load_photoobj(self,tosave=True):
			"""
			If self.dir_obj/PhotoObj.csv exist, load table locally. Otherwise
			download from SDSS and save it locally. 

			Parameters
			-------
			tosave=False: bool
				whether to save the table
			"""
			import  astropy.io.ascii
			# define filename
			filename=self.outer.dir_obj+'PhotoObj.csv'

			if os.path.isfile(filename): 
				print "[obsobj] reading photoobj locally"
				self.photoobj=Table.read(filename,format='ascii.csv',comment='#')
			else:
				print "[obsobj] querying photoobj from SDSS"
				# load photoobj table and store it in obj.sdss.photoobj
				sql="SELECT p.* FROM PhotoObj AS p WHERE p.objid="+str(self.objid)
				#txt=sqlcl.query(sqlc).readlines()
				txt=sqlcl.query(sql).read()
				tabphotoobj=astropy.io.ascii.read(txt)
				print tabphotoobj
				self.photoobj=tabphotoobj
				if tosave and os.path.isdir(os.path.dirname(filename)):
					tabphotoobj.write(filename,format='ascii.csv')




