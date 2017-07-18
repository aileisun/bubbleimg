# hscobj.py
# ALS 2017/05/11

"""
define class hscObj, which can check whether there is such an object and load xid 
"""
import numpy as np
import sys 
import os
from astropy.coordinates import SkyCoord
import astropy.table as at
from astropy.io import fits
import astropy.units as u

from hscsspquery import hscSspQuery
from ..plainobj import plainObj

class hscObj(plainObj):
	def __init__(self, **kwargs):
		"""
		load sdss.xid write files 'hsc_xid.csv' automatically.
		if successful, sdss.status == True

		Instruction
		-----------
		One has to set HSC_SSP STARs account username and password as environmental variable to access hsc data
		  $ export HSC_SSP_CAS_USERNAME
		  $ read -s HSC_SSP_CAS_USERNAME
		  $ export HSC_SSP_CAS_PASSWORD
		  $ read -s HSC_SSP_CAS_PASSWORD


		Params
		------
		ra (float)
		dec (float)
		/either
			dir_obj (string)
		/or 
			dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'
		
		rerun = 's16a_wide' (string): which data base to search for
		release_version = 'dr1' (string): which data base to search for
		search_radius = 2.* u.arcsec


		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		rerun (string): e.g., 's16a_wide'
		release_version (string): e.g., 'dr1'
		search_radius (angle quantity object)

		status (whether the xid and photoboj query were successful)

		optional attr (if querry successful):
			xid
			photoobj
		
		"""
		super(self.__class__, self).__init__(**kwargs)
		self.rerun = kwargs.pop('rerun', 's16a_wide')
		self.release_version = kwargs.pop('release_version', 'dr1')
		self.search_radius = kwargs.pop('search_radius', 2.*u.arcsec)

		self.fp_xid = self.dir_obj+'hsc_xid.csv'
		self.status = self.load_xid()


	def load_xid(self):
		"""
		load xid either locally or remotely and add it as attribute self.xid. If query remotely then save the file to self.fp_xid. 

		Params
		------
		self 

		Return
		------
		status (bool): if true then the loading was successful, false if not
		"""
		xid = self._get_xid()

		if xid is not None:
			self.xid = xid
			for col in xid.colnames: 
				# setting xid attributes 
				setattr(self, col, xid[col][0])
			return True
		else: 
			return False


	def _get_xid(self):
		"""
		return xid.
		Read xid locally if self.dir_obj+'hsc_xid.csv' exist. Otherwise query. 

		Parameters
		------
		self: obj
			contains: 
			self.dir_obj, self.ra, self.dec

		Returns:
		------
		xid: table, or False if failed

		Write Output: (optional)
		------
		self.dir_obj+'xid.csv'
		"""

		fn = self.fp_xid

		if os.path.isfile(fn): # retrieve xid locally
			print "[hscObj] reading xid locally"
		else: # download xid from sdss
			print "[hscObj] querying xid from server"
			self.make_dir_obj()	
			sql = _get_sql(ra=self.ra, dec=self.dec, rerun=self.rerun, radius=float(self.search_radius/u.arcsec))
			hscSspQuery(sql, filename_out=fn, release_version=self.release_version)

		if os.path.isfile(fn): # retrieve xid locally
			if os.stat(fn).st_size > 0:
				xid = at.Table.read(fn, format='ascii.csv', comment='#')

				if len(xid)>=1: 
					if len(xid)>1:
						xid = self._resolve_multiple_sources(xid)

					# sanity check 1
					diffra = (round(xid['ra'], 2) != round(self.ra, 2))
					diffdec = (round(xid['dec'], 2) != round(self.dec, 2))
					if diffra or diffdec:
						raise ValueError("local hsc_xid inconsistent with object")
					# sanity check 2
					if (len(xid) != 1) or (xid['detect_is_primary'][0] !='t'):
						raise Exception("[hscObj] something wrong with xid") 

					return xid

				elif len(xid)<1:
					print "[hscObj] no object found"				
					os.remove(fn)
					return None

			else: 
				print "[hscObj] no object found"
				os.remove(fn)
				return None
		else: 
			print "[hscObj] query failed"
			return None


	def _resolve_multiple_sources(self, xid):
		""" return the xid with only the row that is closest to self.ra, dec"""
		print "[hscObj] multiple primary objects found, choose the closest one"
		c = SkyCoord(self.ra, self.dec, 'icrs', unit='deg')
		crows = [SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg') for row in xid]
		
		a = np.array([c.separation(crow).value for crow in crows])
		xid = at.Table(xid[np.argmin(a)])
		if len(xid) > 1:
			raise Exception("[hscobj] more than one closest object is selected")

		return xid



def _get_sql(ra, dec, radius=2, rerun='s16a_wide'):
	"""
	construct sql query 

	Params
	------
	ra (float): ra in deg decimal
	dec (float): dec in deg decimal
	radius (float): search radius in arcsec
	rerun (string): which rerun to use

	"""

	path=os.path.dirname(sys.modules[__name__].__file__)
	if path == '': 
		path ='.'
	localpath = path+'/'

	fn = localpath+'template.sql'

	with open(fn, 'r') as f:
		sql_template=f.read()
	sql = sql_template.format(rerun=rerun, ra=str(ra), dec=str(dec), radius=str(radius))

	return sql

