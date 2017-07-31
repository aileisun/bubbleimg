# sdssobj.py
# ALS 2017/05/11

"""
define class sdssobj, which can check whether there is such an object and load xid 
"""

import os
from astropy.coordinates import SkyCoord
import astroquery
import astropy.table as at
from astropy.io import fits
import astropy.units as u

from ..plainobj import plainObj

class sdssObj(plainObj):
	def __init__(self, **kwargs):
		"""
		load sdss.xid and sdss.photoboj and write files 'sdss_xid.csv', 'photoobj.csv' automatically.
		if successful, sdss.status == True

		Params
		------
		ra (float)

		dec (float)

		/either
			dir_obj (string)
		/or 
			dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'

		data_release = 12:
			which sdss data release to use

		search_radius = 2.* u.arcsec

		toload_photoobj=True (bool): whether to download photoobj

		writefile=True (bool): whether to write xid ( and photoobj if toload_photoobj=True)


		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		data_release (int)
		search_radius (angle quantity object)
		sdssname (string)
		status (whether the xid and photoboj query were successful)

		optional attr (if querry successful):
			xid
			photoboj
		
		"""
		super(self.__class__, self).__init__(**kwargs)
		self.data_release = kwargs.pop('data_release', 12)
		self.search_radius = kwargs.pop('search_radius', 2.*u.arcsec)

		writefile = kwargs.pop('writefile', True)
		toload_photoobj = kwargs.pop('toload_photoobj', True)

		statusxid = self.load_xid(writefile=writefile)

		if toload_photoobj:
			statusphotoobj = self.load_photoobj(writefile=writefile)

		self.status = statusxid


	def load_xid(self, writefile=True):
		"""
		load xid either locally or remotely and add it as attribute self.xid

		Params
		------
		self 
		writefile=True: 
			if true then write loaded xid to file self.dir_obj/'xid.csv', if it does not already exists

		Return
		------
		status (bool): if true then the loading was successful, false if not
		"""
		xid = self._get_xid(writefile=writefile)

		if xid is not None:
			self.xid = xid
			for col in xid.colnames: 
				# setting xid attributes 
				# ra,dec,objid,run,rerun,camcol,field,z,plate,mjd,fiberID,specobjid,run2d,instrument
				setattr(self, col, xid[col][0])
			return True
		else: 
			return False


	def _get_xid(self, writefile=True):
		"""
		return xid.
		Read xid locally if self.dir_obj+'xid.csv' exist. Otherwise query sdss. 

		Parameters
		------
		self: obj
			contains: 
			self.dir_obj, self.ra, self.dec
		writefile=True

		Returns:
		------
		xid: table

		Write Output: (optional)
		------
		self.dir_obj+'xid.csv'
		"""
		# set up
		photoobj_defs = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field']
		specobj_defs=['z', 'plate', 'mjd', 'fiberID', 'specobjid', 'run2d', 'instrument','sciencePrimary']

		# define filename
		filename = self.dir_obj+'sdss_xid.csv'

		if os.path.isfile(filename): # retrieve xid locally
			print "[sdssobj] reading xid locally"
			xid = at.Table.read(filename,format='ascii.csv',comment='#')

			# sanity check
			diffra = (round(xid['ra'], 2) != round(self.ra, 2))
			diffdec = (round(xid['dec'], 2) != round(self.dec, 2))
			if diffra or diffdec:
				raise ValueError("local sdss_xid inconsistent with object")

			return xid

		else: # download xid from sdss
			print "[sdssobj] querying xid from SDSS"
			c = SkyCoord(self.ra, self.dec, 'icrs', unit='deg')
			result = astroquery.sdss.SDSS.query_region(c, spectro=True, photoobj_fields=photoobj_defs, specobj_fields=specobj_defs, data_release=self.data_release, radius=self.search_radius)

			# Retrieving  sdss ids of the spec sciencePrimary 
			if result is not None:
				xid = result[result['sciencePrimary'] == 1]
				if len(xid) == 1:
					print "[sdssobj] science primary object found"
				elif len(xid) > 1:
					print "[sdssobj] multiple science primary object found, choose the closest one"
					cspecs = [SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg') for row in xid]
					print "[sdssobj] science primary object found"
					a = np.array([c.separation(cspec).value for cspec in cspecs])
					xid = at.Table(xid[np.argmin(a)])
				else:
					print "[sdssobj] no science primary found or duplicate"
					xid = None
			else:
				print "[sdssobj] no object found"
				xid = None

			# write xid
			if (xid is not None) and writefile:
				self.make_dir_obj()
				xid.write(filename,format='ascii.csv',comment='#')
			return xid


	def load_photoobj(self, writefile=True):
		"""
		If self.dir_obj/PhotoObj.csv exist, load table locally. Otherwise
		download from SDSS and save it locally. 

		Parameters
		-------
		writefile=False: bool
			whether to save the table

		Return 
		------
		status: true if photoobj successfully loaded as self.photoobj
		"""
		import astropy.io.ascii
		# define filename
		filename = self.dir_obj+'sdss_photoobj.csv'

		if not hasattr(self, 'xid'):
			xidstatus = self.load_xid(writefile=writefile)
		else: 
			xidstatus = True

		if xidstatus:
			if os.path.isfile(filename): 
				print "[sdssobj] reading photoobj locally"
				self.photoobj = at.Table.read(filename, format='ascii.csv',comment='#')
			else:
				print "[sdssobj] querying photoobj from SDSS"
				# load photoobj table and store it in obj.sdss.photoobj
				sql_query = "SELECT p.* FROM PhotoObj AS p WHERE p.objid="+str(self.objid)
				tabphotoobj = astroquery.sdss.SDSS.query_sql(sql_query, data_release=self.data_release)
				self.photoobj = tabphotoobj

				if writefile:
					self.make_dir_obj()
					tabphotoobj.write(filename, format='ascii.csv',comment='#')

			if hasattr(self, 'photoobj') and isinstance(self.photoobj, at.Table):
				return True
			else: 
				return False
		else: 
			return False



	def get_spec(self):
		"""
		If self.dir_obj/spec.fits exist, load spec locally. Otherwise
		download from SDSS and save it locally. 

		Examples
		------
		spectable = obj.sdss.get_spec(obj)[1].data
		spec, lcoord= spectable['flux'], 10.**spectable['loglam']

		Parameters
		----------
		obj (object of obsobs)

		Returns
		------
		sp (hdulist) or None if it fails

		"""
		fn = self.dir_obj+self.get_spec_filename()

		if not os.path.isfile(fn):
			self.make_spec(overwrite=False)

		if os.path.isfile(fn):
			sp = fits.open(fn)
			return sp
		else: 
			return None


	def make_spec(self, overwrite=False):
		"""
		make spectrum by downloading. If overwrite=False then skip if it exists. 

		Params
		------
		self
		overwrite=True (bool)

		Return
		------
		status (bool)

		Write Outputs
		-------------
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
		self.make_dir_obj()		
		fn = self.dir_obj+self.get_spec_filename()

		if self.status: 
			if not os.path.isfile(fn) or overwrite: 
				print "[sdssObj] download spec"
				sp = astroquery.sdss.SDSS.get_spectra(matches=self.xid, data_release=self.data_release)
				if len(sp) == 1: 
					sp=sp[0]
					self.make_dir_obj()
					sp.writeto(fn, overwrite=overwrite)
					return True
				else: 
					raise ValueError("SDSS spec obj not uniquely identified. ")
					return False
			else: 
				print "[sdssObj] skip download spec as file exists"
				return True
		else: 
			return False


	def get_speclcoord(self, wunit=False):
		"""
		to retrieve SDSS spec and lcoord arrays

		Param
		-----
		wunit=False (bool)  to attach unit or not

		Return
		------
		spec (nparray)
		lcoord (nparray)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_lcoord = u.AA
		"""
		spectable = self.get_spec()[1].data
		spec, lcoord = spectable['flux'], 10.**spectable['loglam']

		if not wunit:	
			return spec, lcoord
		else:	
			u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
			u_lcoord = u.AA
			return spec*u_spec, lcoord*u_lcoord


	def get_spec_filename(self):
		return 'spec.fits'

