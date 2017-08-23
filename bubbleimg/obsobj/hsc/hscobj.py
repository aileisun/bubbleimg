# hscobj.py
# ALS 2017/05/11

"""
define class hscObj, which can check whether there is such an object and load xid 
"""
import numpy as np
import sys 
import os
import astropy.coordinates as ac
import astropy.table as at
from astropy.io import fits
import astropy.io.ascii
import astropy.units as u
import ntpath

import hscsspquery
from ..plainobj import plainObj

fn_table_template_sql = 'table_template.sql'
# fn_photoobj_template_sql = 'photoobj_template.sql'
fn_xid_template_sql = 'xid_template.sql'

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
		
		rerun = 's16a_wide2' (string): which data base to search for

		data_release = 'dr1' (string): which data base to search for

		search_radius = 2.* u.arcsec

		overwrite=False (bool): 
			If true, load xid remotely and rewrite local xid.csv. If false, ready local sdss_xid.csv whenever it eixsts, if not, then load remotely and save to local xid.csv. 


		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		rerun (string): e.g., 's16a_wide2'
		data_release (string): e.g., 'dr1'
		search_radius (angle quantity object)

		status: 
			whether the xid query was successful

		optional attr (if querry successful):
			xid
		
		"""
		super(self.__class__, self).__init__(**kwargs)
		self.rerun = kwargs.pop('rerun', 's16a_wide2')
		self.data_release = kwargs.pop('data_release', 'dr1')
		self.search_radius = kwargs.pop('search_radius', 2.*u.arcsec)

		self.fp_xid = self.dir_obj+'hsc_xid.csv'
		self.fp_xid_temp = self.dir_obj+'hsc_xid_temp.csv'
		self.fp_photoobj = self.dir_obj+'hsc_photoobj.csv'

		overwrite = kwargs.pop('overwrite', False)
		self.load_xid(overwrite=overwrite)


	def get_fp_table(self, tab_name):
		return self.dir_obj + 'hsc_{}.csv'.format(tab_name)


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
				# setting xid attributes 
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

		Params
		------
		self

		Returns:
		------
		xid: table, or None if failed

		Write Output: (optional)
		------
		self.dir_obj+'xid.csv'
		"""
		fn_temp = self.fp_xid_temp
		fn = self.fp_xid

		if os.path.isfile(fn_temp):
			os.remove(fn_temp)

		if not os.path.isfile(fn) or overwrite: # download xid from sdss
			print "[hscObj] querying xid from server"
			self.make_dir_obj()	
			sql = _get_xid_sql(ra=self.ra, dec=self.dec, rerun=self.rerun, search_radius=self.search_radius)
			hscsspquery.hscSspQuery_retry(n_trials=20, sql=sql, filename_out=fn_temp, release_version=self.data_release)
			xid = self._read_xid_from_fp_xid_temp_and_delete_file(fn_temp)

			if xid is not None:
				xid.write(fn, format='ascii.csv', overwrite=overwrite)

		else: # retrieve xid locally
			print "[hscObj] reading xid locally"
			xid = at.Table.read(fn, format='ascii.csv', comment='#')
			self._xid_sanity_check(xid)

		return xid

	def _read_xid_from_fp_xid_temp_and_delete_file(self, fn_temp):
		""" read in xid from file and return it, delete fn_temp after reading """
		if os.path.isfile(fn_temp): 
			if os.stat(fn_temp).st_size > 0:
				xid = at.Table.read(fn_temp, format='ascii.csv', comment='#')
				os.remove(fn_temp)
				xid = self._xid_pick_only_closest(xid)
			else: 
				print "[hscObj] no object found - xid file empty"
				os.remove(fn_temp)
				xid = None
		else: 
			print "[hscObj] query failed"
			xid = None
		return xid


	def _xid_pick_only_closest(self, xid):
		""" takes in xid as arg, picks one, and return """
		if len(xid) == 1:
			self._xid_sanity_check(xid)

		elif len(xid) > 1: 
			print "[hscObj] multiple objects found, picking the closest"
			xid = self._resolve_multiple_sources(xid)
			self._xid_sanity_check(xid)

		elif len(xid) < 1:
			print "[hscObj] no object found"
			xid = None

		return xid


	def _xid_sanity_check(self, xid):
		"""
		check whether:
			1. xid has one unique row
			2. ra, dec consistent with self.ra, self.dec (within self.search_radius). 
			3. detect_primary, detect_is_patch_inner, detect_is_tract_inner are all true

		Raised Exception when any of the checks did not pass. Return None. 
		"""
		# 1) 
		if (len(xid) != 1):
			raise Exception("[hscObj] xid not unique")

		# 2)
		cself = ac.SkyCoord(self.ra, self.dec, 'icrs', unit='deg')
		cxid = ac.SkyCoord(xid['ra'], xid['dec'], 'icrs', unit='deg')
		sep = cself.separation(cxid)
		if sep > self.search_radius:
			raise Exception("[hscObj] xid coordinate inconsistent with object")

		# 3) 
		detect_check = [xid[col][0] == 't' for col in ['detect_is_patch_inner', 'detect_is_tract_inner', 'detect_is_primary']]
		if not all(detect_check):
			raise Exception("[hscObj] detect is not primary or inner patch/tract") 


	def _resolve_multiple_sources(self, xid):
		""" return the xid with only the row that is closest to self.ra, dec"""
		print "[hscObj] multiple primary objects found, choose the closest one"
		c = ac.SkyCoord(self.ra, self.dec, 'icrs', unit='deg')
		crows = [ac.SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg') for row in xid]
		
		a = np.array([c.separation(crow).value for crow in crows])
		xid = at.Table(xid[np.argmin(a)])
		return xid


	def load_photoobj(self, band_columns=[], bands=[], tab_name='forced', all_columns=False, overwrite=False):
		"""
		load hsc photometry table either locally or remotely, depending on whether local file exists and if overwrite=True, and add it as attribute self.photoobj. 

		Params
		------

		self 
		band_columns: Either contains the band_columns when called upon or has a default choice of band_columns if found NULL
			Enter the required band_columns:
			Refer Schema Browser -> "https://hscdata.mtk.nao.ac.jp/schema_browser2/"
			Note: STARs account required
		bands: Either contains the bands req. when called upon or has a default choice of bands if found NULL
			Available bands: g, r, i, z and y
		tab_name='forced': which tab_name to load from remote hsc database
		all_columns=True: Generates SQL code such that all the fields from the table are included
		overwrite=False (bool)

		Return
		------
		status (bool): if true then the loading was successful, false if not
		
		"""
		if not hasattr(self, 'xid'):
			self.load_xid()

		if self.status: 
			status = self._download_photoobj(band_columns=band_columns, bands=bands, tab_name=tab_name, all_columns=all_columns, overwrite=overwrite)

			if status:
				photoobj = at.Table.read(self.fp_photoobj, format='ascii.csv', comment='#')
				self.photoobj = photoobj
				return True
			else:
				return False

		else: 
			print("[hscobj] skip loading photoobj as xid is not successfully created")
			return False


	def _download_photoobj(self, band_columns=[], bands=[], all_columns=False, tab_name='forced', overwrite=False):
		"""
		return photoobj.
		Read photoobj locally if self.dir_obj+'hsc_xid.csv' exist. Otherwise query. 

		Params
		------
		band_columns=[]
		bands=[]
		all_columns=False
		tab_name='forced'
		overwrite=False (bool)

		Returns:
		------
		status

		Write Output: (optional)
		------
		self.dir_obj+'xid.csv'
		"""
		fp = self.fp_photoobj

		object_id = self.xid['object_id'][0]
		sql = _get_photoobj_sql(object_id=object_id, band_columns=band_columns, bands=bands, all_columns=all_columns, rerun=self.rerun, tab_name=tab_name)

		header = {'tab_name': tab_name}
		status = self._download_table_kernel(sql, fp, header=header, overwrite=overwrite)

		return status


	def download_table(self, tab_name='photoz_demp', columns=[], overwrite=False):
		""" 
		download table for the obj given 'object_id' from the hsc catalog to file (hsc_{tab_name}.csv) and return status. 

		Params
		------
		tab_name = 'photoz_demp'
		columns=[]
			by default downloading all columns if nothing is specified. 
		overwrite=False

		Return
		------
		status (bool)
		"""
		fp = self.get_fp_table(tab_name=tab_name)

		object_id = self.xid['object_id'][0]
		sql = _get_table_sql(object_id=object_id, tab_name=tab_name, columns=columns, rerun=self.rerun)

		status = self._download_table_kernel(sql, fp, overwrite=overwrite)

		return status


	def _download_table_kernel(self, sql, fp, header={}, overwrite=False):
		""" 
		the working kernel to download hsc table

		Params
		------
		self
		sql (str)
			sql code to query
		fp (str)
			full file path to be saved to
		header = {} (dictionary)
			set to {'header_name': 'header_content'}
			to be added to the first columns of the downloaded table
		overwrite=False (bool)

		Return
		------
		status
		"""
		fn = ntpath.basename(fp)

		if not os.path.isfile(fp) or overwrite:
			print("[hscobj] querying table {} from HSC".format(fn))

			hscsspquery.hscSspQuery_retry(n_trials=20, sql=sql, filename_out=fp, release_version=self.data_release)

			if os.path.isfile(fp) and (os.stat(fp).st_size > 0):
				print("[hscobj] successful")

				if len(header) > 0:
					_add_header_columns_to_table(fp, header)

				status = True

			else: 
				print("[hscobj] querying table {} from HSC failed".format(fn))
				if os.path.isfile(fp):
					os.remove(fp)
				status = False

		else:
			print("[hscobj] skip querying table {} from HSC as file exists".format(fn))
			status = True

		return status


def _get_xid_sql(ra, dec, search_radius=2 * u.arcsec, rerun='s16a_wide2'):
	"""
	construct sql query 

	Params
	------
	ra (float): ra in deg decimal
	dec (float): dec in deg decimal
	search_radius (float): search radius in arcsec
	rerun (string): which rerun to use

	"""
	localpath = _get_local_path()

	fn = localpath+fn_xid_template_sql

	radius_str = str(float(np.array(search_radius.to(u.arcsec))))

	with open(fn, 'r') as f:
		sql_template=f.read()
	sql = sql_template.format(rerun=rerun, ra=str(ra), dec=str(dec), radius=radius_str)

	return sql


def _add_header_columns_to_table(fp, header):
	table = at.Table.read(fp, format='ascii.csv', comment='#')
	for col in header:
		table = at.hstack([at.Table([[header[col]]], names=[col]), table])
	table.write(fp, format='ascii.csv', overwrite=True)


def _get_local_path():
	path = os.path.dirname(sys.modules[__name__].__file__)
	if path == '': 
		path = '.'
	localpath = path+'/'
	return localpath


def _get_table_sql(object_id, tab_name, columns=[], rerun='s16a_wide2', save_sql=False, fn_sql='hsc_sql.txt', ):
	"""
	return sql code to query a table {tab_name} given object_id

	Params
	------
	object_id (int)
	tab_name (str)
	columns = [] (list)
		by default downloading all columns if nothing is specified. 
	rerun = 's16a_wide2'
	save_sql = False (bool)
	fn_sql = 'hsc_sql.txt' (str)

	Return
	------
	sql (str)
	"""
	localpath = _get_local_path()

	fn = localpath + fn_table_template_sql
	sql_columns = _get_table_sql_columns(columns=columns)

	with open(fn, 'r') as f:
		sql_template = f.read()
	sql = sql_template.format(object_id=object_id, tab_name=tab_name, sql_columns=sql_columns, rerun=rerun, )

	if save_sql:
		with open(fn_sql, "w") as text_file:
			text_file.write(sql)

	return sql


def _get_table_sql_columns(columns=[]):
	"""
	return a string of a list of columns to be inserted into the sql code
	
	Params
	------
	columns=[] (list)

	Return
	------
	sql_columns (str)
	"""

	if len(columns) == 0:
		sql_columns = '*'

	else: 
		sql_columns = ",".join(columns)

	return sql_columns


def _get_photoobj_sql(object_id, band_columns=[], bands=[], all_columns=False, rerun='s16a_wide2', tab_name='forced', save_sql=False, fn_sql='hsc_sql.txt'):
	"""
	construct sql query 

	Params
	------
	object_id
	band_columns=[]
	bands=[]
	all_columns=False
	rerun='s16a_wide2'
	save_sql=False
		if True, save the output sql script to fn_sql
	fn_sql='hsc_sql.txt'

	Return
	------
	sql (str)
	"""

	localpath = _get_local_path()

	fn = localpath+fn_table_template_sql
	sql_columns = _get_photoobj_sql_columns(band_columns=band_columns, bands=bands, all_columns=all_columns)

	with open(fn, 'r') as f:
		sql_template=f.read()
	sql = sql_template.format(object_id=object_id, tab_name=tab_name, sql_columns=sql_columns, rerun=rerun, )
	
	if save_sql:
		with open(fn_sql, "w") as text_file:
			text_file.write(sql)

	return sql


def _get_photoobj_sql_columns(band_columns=[], bands=[], all_columns=False):
	
	"""
	Creates SQL commands for required band_columns in the form of a single string array
	
	Params
	------
	band_columns=[]
	bands=[]
	all_columns=False


	Return
	------
	Returns a string containing the SQL commands for data acquisition
	"""

	default = 'object_id, ra, dec, patch_id, tract, patch, patch_s, parent_id, deblend_nchild, detect_is_patch_inner, detect_is_tract_inner, detect_is_primary, '

	if all_columns:
		sql_columns = '*'

	else:
		sql_columns = ""
		if len(band_columns) == 0:
			band_columns = ['mag_kron', 'mag_kron_err', 'flux_kron_flags', 'flux_kron_radius', 'mag_aperture10', 'mag_aperture15']

		if len(bands) == 0:
			bands = ['g', 'r', 'i', 'z', 'y'] 

		count = 0
		for x in bands:
			for y in band_columns:
				count+=1
				if count!=len(bands)*len(band_columns):
					sql_columns += x + y + ","
				else:
					sql_columns += x + y
		sql_columns = default + sql_columns

	return sql_columns

