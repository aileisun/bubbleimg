# batch.py
# ALS 2017/05/29

import numpy as np
import astropy.table as at
from astropy.io import ascii
import os
import sys
import shutil
import copy
import multiprocessing as mtp

from .. import obsobj
import mtp_tools

class Batch(object):
	def __init__(self, survey, obj_naming_sys='sdss', args_to_list=[], **kwargs):
		"""
		Batch

		dir_batch is not automatically created, please call mkdir_batch()

		Params
		----------

		/either
			dir_batch (string): path of the directory
		/or 
			dir_parent (string): attr dir_batch is set to dir_parent+name+'/'
			name (string): name of the batch

		/either
			catalog (astropy table): 
				list of objects with columns ['ra', 'dec', ...] or ['RA', 'DEC', ...]
		/or 
			fn_cat (str):
				path to the catalog table file
		/or
			nothing: if dir_batch/list.csv exist, then read in list.csv as catalog 

		survey = (str): either 'sdss' or 'hsc'
		obj_naming_sys = 'sdss' (str): how to name objects in the batch
		args_to_list = [] (list): 
			list of names of additional arguments to pass to list from catalog, e.g., ['z']. 
				

		Attributes
		----------
		dir_batch (str)
		name (str)
		catalog (astropy table): the input catalog
		survey (str)
		obj_naming_sys = 'sdss' (str)

		list (astropy table):
			list of objects with columns ['ra', 'dec', 'obj_name']. 
		list_good (astropy table): 
			same as list but contains only good object
		list_except (astropy table):
			contain objects with exceptions
		args_to_list (list of str):
			list of argument names to include in lists

		dir_good (str)
		dir_except (str)
		fp_list_good (str)
		fp_list_except (str)
		fp_list (str)
		"""

		# set dir_batch, name
		if 'dir_batch' in kwargs:
			self.dir_batch = kwargs.pop('dir_batch', None)
			self.name = self.dir_batch.split('/')[-2]

			if 'name' in kwargs:
				if kwargs['name'] != self.name:
					raise Exception("[batch] input name inconsistent with dir_batch")

		elif ('dir_parent' in kwargs) and ('name' in kwargs):
			self.dir_parent = kwargs.pop('dir_parent', None)
			self.name = kwargs.pop('name', None)
			self.dir_batch = self.dir_parent+self.name+'/'

		else:
			raise Exception("[batch] dir_batch or dir_parent/name not specified")

		# set directories
		self.fp_list = self.dir_batch+'list.csv'
		self.dir_good = self.dir_batch+'good/'
		self.dir_except = self.dir_batch+'except/'
		self.fp_list_good = self.dir_good+'list_good.csv'
		self.fp_list_except = self.dir_except+'list_except.csv'

		# set catalog
		if 'catalog' in kwargs:
			self.catalog = kwargs.pop('catalog', None)

		elif 'fn_cat' in kwargs:
			fn_cat = kwargs.pop('fn_cat', None)
			if os.path.splitext(fn_cat)[1] == '.fits':
				self.catalog = at.Table.read(fn_cat)
			elif os.path.splitext(fn_cat)[1] == '.csv':
				self.catalog = at.Table.read(fn_cat, format='ascii.csv', comment='#')
			else:
				raise Exception("[batch] input catalog file extension not recognized")

		elif os.path.isfile(self.fp_list):
			self.catalog = at.Table.read(self.fp_list)
			args_to_list = [col for col in self.catalog.colnames if col not in ['ra', 'dec', 'obj_name']]

		else:
			raise Exception("[batch] input catalog not specified")

		# setting other attribute
		self.survey = survey
		self.obj_naming_sys = obj_naming_sys
		self.args_to_list = args_to_list
		self._set_attr_list()
		self._set_attr_list_good()
		self._set_attr_list_except()

		# sanity check
		if self.dir_batch[-1] != '/':
			raise Exception("[batch] dir_batch not a directory path")

		if self.survey not in ['sdss', 'hsc']:
			raise Exception("[batch] survey not recognized")


	def mkdir_batch(self):
		if not os.path.isdir(self.dir_batch):
			os.makedirs(self.dir_batch)


	def iterlist(self, func, listargs=[], listname='good', overwrite=False, processes=None, **kwargs): 
		"""
		Apply function to each of the objects in list (default: good). Nothing is done to the function return. 

		Params
		------
		func (function):
		    a function that takes in parameter (obj, overwrite, **kwargs), where kwargs includes listargs, and returns result (usually status [bool]). 

		listargs=[] (list of str): names of arguments to pass to func from self.list

		listname='good' (str): name of the list to run on, either 'good' or 'except'. 

		overwrite=False

		processes = None (int):
			How many processes to use for multiprocessing. Default is None, the default of multiprocessing.Pool. If processes == -1, then it will be ran sequentially and no multiprocessing is used. 

		**kwargs:
	    	other arguments to be passed to func

	    Return
	    ------
	    results (list):
	    	a list of all the results returned by fun, usually "status" (bool array)
		"""
		self._check_folders_consistent_w_list()

		dir_list = self._get_dir_list_from_listname(listname=listname)
		lst = self._get_list_of_listname(listname=listname)

		if processes != -1: # to run multiprocessing
			ikernel_partial = mtp_tools.partialmethod(self._iterlist_kernel, func=func, listargs=listargs, dir_list=dir_list, overwrite=overwrite, **kwargs)
			p = mtp.Pool(processes=processes)
			results = p.map(ikernel_partial, lst)
			p.close()
			p.join()

		else:  # to run sequantially
			results = []
			for row in lst:
				result = self._iterlist_kernel(row=row, func=func, listargs=listargs, dir_list=dir_list, overwrite=overwrite, **kwargs)
				results += [result]

		return results


	def _iterlist_kernel(self, row, func, listargs, dir_list, overwrite, **kwargs):
		""" the kernel to be iterated over (with different input row) in self.iterlist() """
		ra = row['ra']
		dec = row['dec']
		obj_name = row['obj_name']

		print("[batch] {obj_name} iterating".format(obj_name=obj_name))

		for arg in listargs:
			kwargs.update({arg: row[arg]})

		obj = obsobj.obsObj(ra=ra, dec=dec, dir_parent=dir_list, obj_naming_sys=self.obj_naming_sys, overwrite=overwrite)
		obj.survey = self.survey		

		result = func(obj, overwrite=overwrite, **kwargs)
		return result


	def get_ith_obj_from_list(self, iobj, listname='good'):
		""" 
		return the obj for the ith object from the specified list (good or except)
		"""
		dir_list = self._get_dir_list_from_listname(listname=listname)
		lst = self._get_list_of_listname(listname=listname)

		ra = lst['ra'][iobj]
		dec = lst['dec'][iobj]
		obj_name = lst['obj_name'][iobj]

		obj = obsobj.obsObj(ra=ra, dec=dec, dir_parent=dir_list, obj_naming_sys=self.obj_naming_sys)

		return obj


	def compile_table(self, fn_tab, overwrite=False):
		"""
		compile table fn_tab for the entire batch and creates self.dir_batch+fn_tab.

		Params
		------
		self
		fn_tab (str):
			file name of the table to be compiled, e.g., 'sdss_photoobj.csv'.
		overwrite (bool)

		Return
		------
		status (bool)

		Write output
		------------
		self.dir_batch + fn_tab:
			e.g., 'batch_rz/sdss_photoobj.csv'
		"""		
		self._check_folders_consistent_w_list()

		fp = self.dir_batch+fn_tab

		if not os.path.isfile(fp) or overwrite:
			print("[batch] compiling table {}".format(fn_tab))

			if len(self.list_good) > 0:
				obj1 = self.get_ith_obj_from_list(iobj=0, listname='good')
				header = _extract_line_from_file(obj1.dir_obj+fn_tab, iline=0)

				lines_data = self.iterlist(func=self._iterfunc_extract_line_from_file, listargs=[], listname='good', overwrite=False, **{'fn': fn_tab})
				tab_data = ascii.read([header]+lines_data)

				self._rename_list_args(tab_data)
				tab_good = at.hstack([self.list_good, tab_data])
			else:
				tab_good = at.Table()

			if len(self.list_good) > 0 and len(self.list_except) > 0:
				row_masked = at.Table(tab_data[0], masked=True)
				row_masked.mask = True
				tab_masked = at.vstack([row_masked for i in range(len(self.list_except))])
				self._rename_list_args(tab_masked)
				tab_except = at.hstack([self.list_except, tab_masked])
			else:
				tab_except = at.Table()

			tab = at.vstack([tab_good, tab_except])

			if len(tab) > 0:
				tab.sort('ra')
				tab.write(fp, overwrite=overwrite)
			else:
				# in case if there is no good object
				print("[batch] skipped compiling table {} as no data to compile".format(fn_tab))

		else:
			print("[batch] skipped compiling table {} as file exists".format(fn_tab))

		status = os.path.isfile(fp)	
		return status 


	def steal_columns(self, tab, colnames=['z'], keys=['ra', 'dec']):
		"""
		steal columns (colnames) from a table (tab) and add them to list (including list_good, list_except). The table has to have identical objects and keys (ra, dec) as the list in batch. 

		One needs to run compile_table to update compiled tables. 

		Params
		------
		self
		tab (astropy table):
			a table with identical objects and keys (ra, dec) as the list
		colnames=['z']:
			a list of columns names to steal (string)
		keys=['ra', 'dec']:
			a list of the keys used for matching. 

		Write output
		------------
		rewrites list, list_good, and list_except, and update then as self.*. 
		"""

		for col in colnames:
			if col in self.list.colnames:
				raise Exception("[batch] column to steal already exists in list")

		cols_trim = keys + colnames
		tab_trim = tab[cols_trim]

		for lst in [self.list, self.list_good, self.list_except]:
			if len(lst) > 0:
				lst_joined = at.join(lst, tab_trim, keys=keys, join_type='left')
				columns_toadd = [lst_joined[cn] for cn in colnames]
				lst.add_columns(columns_toadd)
			else:
				columns_toadd = [tab_trim[cn] for cn in colnames]
				lst.add_columns([at.Column(name = c.name, dtype=c.dtype, meta=c.meta) for c in columns_toadd])

		self._write_all_lists()


	def remove_columns(self, colnames=['z']):
		"""
		remove columns (colnames) from a list, list_good, and list_except

		Params
		------
		self
		colnames=['z']:
			a list of columns names to steal (string)

		Write output
		------------
		rewrites list, list_good, and list_except, and update then as self.*. 
		"""

		for col in colnames:
			if col not in self.list.colnames:
				raise Exception("[batch] column {} to remove not in list".format(col))

		for lst in [self.list, self.list_good, self.list_except]:
			lst.remove_columns(colnames)

		self._write_all_lists()


	def _batch__build_core(self, func_build, overwrite=False, processes=None, **kwargs):
		"""
		Build batch by making directories and running func_build. To be called by child class. 
		Good objects will be stored in dir_batch/good/.
		Except objects will be stored in dir_batch/except/.

		Params
		------
		func_build=self._func_build (funcion):
			A funciton that takes (obj, overwrite, **kwargs) as param, 
			and returns status (bool)
			The obj contains: ra, dec, dir_parent, dir_obj, name, overwrite, and "survey" of batch
		overwrite=False (bool)
		processes = None (int):
			How many processes to use for multiprocessing. Default is None, the default of multiprocessing.Pool. If processes == -1, then it will be ran sequentially and no multiprocessing is used. 
		**kwargs:
			 to be entered into func_build() in the kwargs part

		Return
		----------
		status: True if successful. 
		"""
		self.mkdir_batch()

		for directory in [self.dir_good, self.dir_except]:
			if not os.path.isdir(directory):
				os.makedirs(directory)

		if overwrite:
			# reset list_good and list_except
			self.list_good = self._create_empty_list_table()
			self.list_except = self._create_empty_list_table()

		self._write_all_lists()

		if processes != -1: # to run multiprocessing
		
			bkernel_partial = mtp_tools.partialmethod(self._buildcore_kernel, func_build=func_build, overwrite=overwrite, **kwargs)
			p = mtp.Pool(processes=processes)
			results = p.map(bkernel_partial, self.list)
			p.close()
			p.join()

		else:  # to run sequantially
			for row in self.list:
				self._buildcore_kernel(row=row, func_build=func_build, overwrite=overwrite, **kwargs)

		self._set_attr_list_good()
		self._set_attr_list_except()
		self._write_a_list(listname='good')
		self._write_a_list(listname='except')
		self._check_folders_consistent_w_list()

		# status reflects that all is ran (list = list_good + list_except)
		list_ran = []
		for lst in [self.list_good, self.list_except]:
			if len(lst)>0:
				list_ran += list(lst['obj_name'])

		obj_names = self.list['obj_name']
		status = (obj_names.sort() == (list_ran).sort())

		if status: 
			print("[batch] building batch {0} done".format(self.name))
		else: 
			print("[batch] building batch {0} unfinished".format(self.name))

		return status


	def _buildcore_kernel(self, row, func_build, overwrite, **kwargs):
		""" the kernel to be iterated over (with different input row) in self._batch__build_core() """
		ra = row['ra']
		dec = row['dec']
		obj_name = row['obj_name']
		dir_obj_good = self.dir_good+row['obj_name']+'/'
		dir_obj_except = self.dir_except+row['obj_name']+'/'

		if (obj_name not in self.list_good['obj_name']) and (obj_name not in self.list_except['obj_name']):
			print("[batch] {obj_name} building".format(obj_name=obj_name))

			try:
				obj = obsobj.obsObj(ra=ra, dec=dec, dir_parent=self.dir_good, obj_naming_sys=self.obj_naming_sys, overwrite=overwrite)
				obj.survey = self.survey
				status = func_build(obj=obj, overwrite=overwrite, **kwargs)

			except KeyboardInterrupt, e:
				print("[batch] func_build() encounters exception {0}".format(str(e)))
				if os.path.isdir(dir_obj_good) and (obj_name not in self.list_good['obj_name']):
					shutil.rmtree(dir_obj_good)
				sys.exit(1)

			if status:
				print("[batch] Successful")
				with open(self.fp_list_good, "a") as f:
					ascii.write(row, output=f, format='no_header', delimiter=',')

			else: 
				print("[batch] Failed, moving to except/. ")
				if os.path.isdir(dir_obj_except):
					shutil.rmtree(dir_obj_except)
				shutil.move(dir_obj_good, dir_obj_except)
				with open(self.fp_list_except, "a") as f:
					ascii.write(row, output=f, format='no_header', delimiter=',')

		else: 
			print("[batch] {obj_name} skipped".format(obj_name=obj_name))


	def _set_attr_list(self):
		"""
		extract list (table of cols ['ra', 'dec', 'obj_name']) from catalog
		"""
		if ('ra' in self.catalog.colnames) and ('dec' in self.catalog.colnames):
			self.list = self.catalog['ra', 'dec']
		elif ('RA' in self.catalog.colnames) and ('DEC' in self.catalog.colnames):
			self.list = self.catalog['RA', 'DEC']
			self.list.rename_column('RA', 'ra')
			self.list.rename_column('DEC', 'dec')
		else: 
			raise NameError("[batch] catalog table does not contain ra, dec or RA, DEC")

		objname = at.Column(name='obj_name', dtype='S64', length=len(self.list))
		for i, row in enumerate(self.list):
			ra = row['ra']
			dec = row['dec']
			objname[i] = obsobj.objnaming.get_obj_name(ra=ra, dec=dec, obj_naming_sys=self.obj_naming_sys)

		self.list.add_column(objname)

		if len(self.args_to_list) > 0:
			self.args_to_list_dtype = self.catalog[self.args_to_list].dtype
			for arg in self.args_to_list:
				self.list[arg] = self.catalog[arg]
		else: 
			self.args_to_list_dtype = None

		self.list.sort('ra')


	def _set_attr_list_good(self):
		""" read sorted list_good from file and set it to attribute """
		self.list_good = self._read_a_list_sorted(listname='good')


	def _set_attr_list_except(self):
		""" read sorted list_except from file and set it to attribute """
		self.list_except = self._read_a_list_sorted(listname='except')


	def _create_empty_list_table(self):
		"""
		Return an empty table with columns ra, dec, obj_name and those in args_to_list 
		"""
		d = np.dtype([('ra', 'float64'), ('dec', 'float64'), ('obj_name', 'S64')])
		tab0 = at.Table(dtype=d)

		if len(self.args_to_list) > 0:
			tab1 = at.Table(dtype=self.args_to_list_dtype)
			tab =  at.hstack([tab0, tab1])
		else:
			tab = tab0

		return tab


	def _read_a_list_sorted(self, listname=''):
		""" 
		Read and return a list. the list is sorted by ra before returning

		Params
		------
		listname (str):
			'' for list, 'good' for list_good, and 'except' for list_except
		"""
		fn = self._get_fp_of_listname(listname=listname)
		if os.path.isfile(fn):
			lst = at.Table.read(fn)
			if len(lst) == 0:
				lst = self._create_empty_list_table()
		else:
			lst = self._create_empty_list_table()
		lst.sort('ra')
		return lst


	def _write_a_list(self, listname=''):
		""" 
		write a list from attribute to file

		Params
		------
		listname (str):
			'' for list, 'good' for list_good, and 'except' for list_except
		"""
		fn = self._get_fp_of_listname(listname=listname)
		lst = self._get_list_of_listname(listname=listname)
		self.mkdir_batch()
		lst.write(fn, format='ascii.csv', overwrite=True)


	def _write_all_lists(self):
		""" write all the lists from attributes to file """
		for listname in ['', 'good', 'except']:
			self._write_a_list(listname=listname)


	def _get_dir_list_from_listname(self, listname=''):
		if listname == '':
			dir_fp = self.dir_ba
		elif listname == 'good':
			dir_fp = self.dir_good
		elif listname == 'except':
			dir_fp = self.dir_except
		return dir_fp


	def	_get_fp_of_listname(self, listname=''):
		if listname == '':
			fp = self.fp_list
		elif listname == 'good':
			fp = self.fp_list_good
		elif listname == 'except':
			fp = self.fp_list_except

		return fp

	def _get_list_of_listname(self, listname=''):
		if listname == '':
			lst = self.list
		elif listname == 'good':
			lst = self.list_good
		elif listname == 'except':
			lst = self.list_except

		return lst


	def _iterfunc_extract_line_from_file(self, obj, fn, iline=1, overwrite=False):
		""" 
		function to be iterated for each object to return a particular line in the file as a string

		Params
		------
		self
		obj
		fn (str):
			file name 
		iline (int):
			the number of the row
		overwrite (bool):
			can be ignored
		"""
		return _extract_line_from_file(obj.dir_obj+fn, iline=iline)


	def _rename_list_args(self, tab):
		""" rename the arguments that could conflict with the ones in list """

		args = ['ra', 'dec', 'obj_name'] + self.args_to_list

		for arg in args:
			if arg in tab.colnames:
				tab.rename_column(arg, arg+'_1')


	def _check_folders_consistent_w_list(self):
		"""
		check that the created obj directories is the same as the list, for list, list_good, and list_except. 
		"""	
		dp_good = self.dir_good
		dp_excp = self.dir_except

		lfolder_good = [obj_name for obj_name in os.listdir(dp_good) if os.path.isdir(os.path.join(dp_good, obj_name))]
		lfolder_excp = [obj_name for obj_name in os.listdir(dp_excp) if os.path.isdir(os.path.join(dp_excp, obj_name))]
		lfolder = lfolder_good + lfolder_excp

		for name, thelist, thefolders in (('list', self.list, lfolder), ('good', self.list_good, lfolder_good), ('except', self.list_except, lfolder_excp)):

			arr_list = np.array(thelist['obj_name'])
			arr_fold = np.array(thefolders)

			arr_list.sort()
			arr_fold.sort()

			n_list = len(arr_list)
			n_fold = len(arr_fold)

			if any([n_list>0, n_fold>0]):
				if n_list != n_fold:
					raise Exception("[batch] number of object folders {n_fold} inconsistent with the list ({name}, {n_list}) in the batch".format(n_fold=n_fold, name=name, n_list=n_list))

				if not all(arr_list == arr_fold): 
					raise Exception("[batch] list of object folders inconsistent with the list in the batch")


def _extract_line_from_file(fn, iline=1, comment='#'): 
	""" return the iline-th line of the file which is non-empty and does not start with the comment # """

	with open(fn, 'r') as f:
		data = f.read()
	lines = data.split('\n')

	lines_noncomment = []
	for line in lines:
		if len(line) > 0:
			if (line[0] != comment):
				lines_noncomment += [line]

	return lines_noncomment[iline]


