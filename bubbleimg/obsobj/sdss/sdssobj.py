# sdssobj.py
# ALS 2017/05/11

"""
define class sdssobj, which can check whether there is such an object and load xid 
"""

import os
import numpy as np
import astropy.coordinates as ac
import astroquery
import astropy.table as at
from astropy.io import fits
import astropy.units as u
import astropy.io.ascii
import requests

from ..plainobj import plainObj

class sdssObj(plainObj):
	def __init__(self, **kwargs):
		"""
		create sdssobj
		load sdss.xid locally or remotely, depending on whether dir_obj/sdss_xid.csv exists locally or if overwrite == True. 
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

		science_primary=True (bool)
			If true, then only sciencePrimary spectrum is taken

		overwrite=False (bool): 
			If true, load xid remotely and rewrite local xid.csv. If false, ready local sdss_xid.csv whenever it eixsts, if not, then load remotely and save to local xid.csv. 

		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		data_release (int)
		search_radius (angle quantity object)
		science_primary=True (bool)
		sdssname (string)
		status (whether the xid and photoboj query were successful)

		optional attr (if querry successful):
			xid
			photoboj
		
		"""
		super(self.__class__, self).__init__(**kwargs)
		self.data_release = kwargs.pop('data_release', 12)
		self.search_radius = kwargs.pop('search_radius', 2.*u.arcsec)
		self.science_primary = kwargs.pop('science_primary', True)

		self.fp_photoobj = self.dir_obj+'sdss_photoobj.csv'
		self.fn_spec = self.dir_obj+'spec.fits'

		overwrite = kwargs.pop('overwrite', False)
		self.load_xid(overwrite=overwrite)


	def load_xid(self, overwrite=False):
		"""
		load xid either locally or remotely and add it as attribute self.xid

		Params
		------
		self 
		overwrite=False: 
			If true then always load xid remotely and rewrites local file "sdss_xid.csv". Otherwise, read locally whenever possible. 

		Return
		------
		status (bool): if true then the loading was successful, false if not
		"""
		xid = self._get_xid(overwrite=overwrite)

		if xid is not None:
			self.xid = xid
			for col in xid.colnames: 
				setattr(self, col, xid[col][0])

			status = True
		else: 
			status = False

		self.status = status
		return status


	def _get_xid(self, overwrite=False):
		"""
		return xid. 
		If overwrite == true then always load xid remotely and rewrites local file "sdss_xid.csv". Otherwise, read locally whenever file exists or load remotely and write file if not. 

		Parameters
		------
		self: obj
			contains: 
			self.dir_obj, self.ra, self.dec
		overwrite=False

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
		fn = self.dir_obj+'sdss_xid.csv'

		if not os.path.isfile(fn) or overwrite:
			# download xid from sdss
			print "[sdssobj] querying xid from SDSS"
			c = ac.SkyCoord(self.ra, self.dec, 'icrs', unit='deg')

			func_query = astroquery.sdss.SDSS.query_region
			kwargs_query = dict(coordinates=c, spectro=True, photoobj_fields=photoobj_defs, specobj_fields=specobj_defs, data_release=self.data_release, radius=self.search_radius)
			result = _retry_sdss_query(func_query, n_trials=5, **kwargs_query)

			# Retrieving  sdss ids of the spec sciencePrimary 
			if result is not None:
				if self.science_primary:
					xid = result[result['sciencePrimary'] == 1]
				else: 
					xid = result

				if xid is not None:
					if len(xid) == 1:
						print "[sdssobj] science primary object found"
					elif len(xid) > 1:
						print "[sdssobj] multiple science primary object found, choose the closest one"
						cspecs = [ac.SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg') for row in xid]
						print "[sdssobj] science primary object found"
						a = np.array([c.separation(cspec).value for cspec in cspecs])
						xid = at.Table(xid[np.argmin(a)])
					else:
						print "[sdssobj] no science primary found or duplicate"
						xid = None
				else:
					pass
			else:
				print "[sdssobj] no object found"
				xid = None

			# write xid
			if (xid is not None):
				self._xid_sanity_check(xid)
				self.make_dir_obj()
				xid.write(fn, format='ascii.csv', comment='#', overwrite=overwrite)

		else: 
			# retrieve xid locally
			print "[sdssobj] reading xid locally"
			xid = at.Table.read(fn, format='ascii.csv', comment='#')

			if (xid is not None):
				self._xid_sanity_check(xid)

		return xid


	def _xid_sanity_check(self, xid):
		"""
		check whether:
			1. xid has one unique row
			2. ra, dec consistent with self.ra, self.dec (within self.search_radius). 

		Raised Exception when any of the checks did not pass. Return None. 
		"""
		# 1) 
		if (len(xid) != 1):
			raise Exception("[sdssObj] xid not unique")

		# 3)
		cself = ac.SkyCoord(self.ra, self.dec, 'icrs', unit='deg')
		cxid = ac.SkyCoord(xid['ra'], xid['dec'], 'icrs', unit='deg')
		sep = cself.separation(cxid)
		if sep > self.search_radius:
			raise Exception("[sdssObj] xid coordinate inconsistent with object")


	def load_photoobj(self, overwrite=False):
		"""
		load sdss photometry table either locally or remotely, depending on whether local file exists and if overwrite=True, and add it as attribute self.photoobj. 

		Params
		------
		self
		overwrite=False (bool)

		Return 
		------
		status: true if photoobj successfully loaded as self.photoobj
		"""
		# define filename

		if not hasattr(self, 'xid'):
			self.load_xid(overwrite=overwrite)

		if self.status:
			photoobj = self._get_photoobj(overwrite=overwrite)

			if photoobj is not None:
				self.photoobj = photoobj
				return True
			else:
				return False
		else: 
			return False


	def _get_photoobj(self, overwrite=True):
		"""
		return photoobj.
		Read photoobj locally if self.dir_obj+'sdss_xid.csv' exist. Otherwise query. 

		Params
		------
		overwrite=False (bool)

		Returns:
		------
		photobj: table, or None if failed

		Write Output: (optional)
		------
		self.fp_photoobj
		"""

		fn = self.fp_photoobj

		if not os.path.isfile(fn) or overwrite:
			print "[sdssobj] querying photoobj from SDSS"
			sql_query = "SELECT p.* FROM PhotoObj AS p WHERE p.objid="+str(self.objid)
			# photoobj = astroquery.sdss.SDSS.query_sql(sql_query, data_release=self.data_release)

			func_query = astroquery.sdss.SDSS.query_sql
			kwargs_query = dict(sql_query=sql_query, data_release=self.data_release)
			photoobj = _retry_sdss_query(func_query, n_trials=5, **kwargs_query)

			if photoobj is not None:
				if len(photoobj) > 0 :
					self.make_dir_obj()
					photoobj.write(fn, format='ascii.csv', comment='#', overwrite=overwrite)
				else:
					photoobj = None

		else: 
			print "[sdssobj] reading photoobj locally"
			photoobj = at.Table.read(fn, format='ascii.csv',comment='#')

		return photoobj


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
		fn = self.fn_spec

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
		fn = self.fn_spec

		if self.status: 
			if not os.path.isfile(fn) or overwrite: 
				print "[sdssObj] download spec"
				# sp = astroquery.sdss.SDSS.get_spectra(matches=self.xid, data_release=self.data_release)
				func_query = astroquery.sdss.SDSS.get_spectra
				kwargs_query = dict(matches=self.xid, data_release=self.data_release)
				sp = _retry_sdss_query(func_query, n_trials=5, **kwargs_query)

				if sp is not None:
					if len(sp) == 1: 
						sp=sp[0]
						self.make_dir_obj()
						sp.writeto(fn, overwrite=overwrite)
						status = True
					else: 
						raise ValueError("SDSS spec obj not uniquely identified. ")
						status = False
				else:
					status = False
			else: 
				print "[sdssObj] skip download spec as file exists"
				status = True
		else: 
			status = False

		return status


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


def _retry_sdss_query(func_query, n_trials=5, **kwargs_query):
	"""
	run sdss query and retries for up to n_trials times if requests exceptions are raised, such as ConnectionErrors. Input the desired request using func_query and **kwargs_query. 
	"""

	for _ in range(n_trials):
		try:
			results = func_query(**kwargs_query)
			return results
			break
		except Exception as e:
			print("[sdssobj] retrying as error detected: "+str(e))


		# except requests.exceptions.RequestException as e:
		# 	print("[sdssobj] retrying as error detected: "+str(e))


