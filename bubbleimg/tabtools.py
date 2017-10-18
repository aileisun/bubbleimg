import io
import os
import numpy as np
import astropy.table as at
from astropy.io import ascii

def write_row(fn, row, condi, overwrite=False, append=False):
	"""
	write row (append) to file. If the row already exists, according to the condi conditions, then this row is overwritten (or not) depending on the overwrite parameter. If append = True then write rows to the end without deleting previous duplications. 

	Params
	------
	fn (str)
	row (astropy tab)
	condi (dictionary)
		e.g., condi = {'imgtag': 'OIII5008_I'}
	overwrite=False
	append=False
	"""
	withheader = not os.path.isfile(fn)

	if append:
		append_row_to_end(fn, row, withheader=withheader)

	elif overwrite:
		fn_delete_row(fn, condi)
		append_row_to_end(fn, row, withheader=withheader)

	elif (not fn_has_row(fn, condi)):
		append_row_to_end(fn, row, withheader=withheader)

	else: 
		print("[tabtools] skip writing row as it exists")


def append_row_to_end(fn, row, withheader=False):
	""" append the row to the end of file """
	rowstring = tab_to_string(row, withheader=withheader)

	with open(fn, 'a') as f_to:
		f_to.write(rowstring)


def tab_to_string(tab, withheader=False):
	""" turn table into string with each line seperated by \n """
	with io.BytesIO() as f_temp: 
		tab.write(f_temp, format='ascii.csv')
		tabstring = f_temp.getvalue()

	if not withheader:
		# take out header
		tabstring = '\n'.join(tabstring.splitlines()[1:]) + '\n'

	return tabstring


def fn_has_row(fn, condi):
	""" 
	return if table has a line with column (key) equals to value. If file does not exist return false. 

	Params
	------
	fn (str): file name
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

		e.g., condi = {'imgtag': 'OIII5008_I'}
	"""
	if os.path.isfile(fn):
		tab = at.Table.read(fn)
		result = tab_has_row(tab, condi)

	else: 
		result = False

	return result


def tab_has_row(tab, condi):
	""" 
	return if table has a line with column (key) equals to value. 

	Params
	------
	tab: table
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

		e.g., condi = {'imgtag': 'OIII5008_I'}
	"""
	select = get_select(tab, condi)
	return np.sum(select) > 0


def fn_delete_row(fn, condi):
	""" 
	delete the lines in table that satisfies the condition
	Params
	------
	fn (str): file name
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

		e.g., condi = {'imgtag': 'OIII5008_I'}
	"""
	if os.path.isfile(fn):
		tab = at.Table.read(fn)
		tab = tab_delete_row(tab, condi)
		tab.write(fn, overwrite=True)


def tab_delete_row(tab, condi):
	""" 
	delete the lines in table that satisfies the condition

	Params
	------
	tab: table
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

		e.g., condi = {'imgtag': 'OIII5008_I'}

	Return
	------
	tab
	"""
	select = get_select(tab, condi)

	if np.sum(select)>0:
		tab.remove_rows(select)

	return tab


def tab_extract_row(tab, condi):
	"""
	return a table of only the extracted rows that meet the condition.

	Params
	------
	tab: table
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

	Return
	------
	tab
	"""
	select = get_select(tab, condi)
	return tab[select]

	
def get_select(tab, condi):
	""" return boolean array indicating whether each row of the tab is selected """
	select_arr = [[str(tab[key][i]) == str(condi[key]) for i in range(len(tab))] for key in condi]
	select = np.all(select_arr, axis=0)
	return select


def summarize(fn_in, fn_out, columns=[], condi={}, overwrite=False):
	"""
	Summarize the table 'fn_in' and write the results to 'fn_out'. 
	For each of the column in columns, take the mean, std, median, and 16%, 84% quantile. 
	All the other columns that are not specified in columns and condi are ignored. 

	Params
	------
	fn_in
	fn_out
	columns=[] (list of string)
		list of column names, e.g., ['area_ars', 'dmax_ars']. Default: all columns. 
	condi={} 
		conditions, e.g. {'imgtag': 'OIII5008_I'}
	overwrite=False

	Return
	------
	status (bool)
	"""
	if not os.path.isfile(fn_out) or overwrite:
		tab_in = at.Table.read(fn_in)

		if len(condi)>0:
			tab_select = tab_extract_row(tab_in, condi=condi)
			tab_sum = tab_select[condi.keys()][0] # creating headers
		else: 
			tab_select = tab_in
			tab_sum = at.Table() # no headers

		if len(columns)==0:
			columns = tab_in.colnames

		# calculation
		for col in columns:
			if not col in condi.keys():
				arr = tab_select[col]
				if arr.dtype in [float, int]:
					var_mean = np.mean(arr)
					var_std = np.std(arr)
					var_median = np.median(arr)
					var_p16 = np.percentile(arr, 16)
					var_p84 = np.percentile(arr, 84)

					tab_stat = at.Table([[var_mean], [var_std], [var_median], [var_p16], [var_p84], ], names=[col+tag for tag in ['_mean', '_std', '_median', '_p16', '_p84']])
					tab_sum = at.hstack([tab_sum, tab_stat])

		tab_sum.write(fn_out, overwrite=overwrite)

	else: 
		print("[tabtools] skip summarizing as files exist")

	return os.path.isfile(fn_out)


def extract_line_from_file(fn, iline=1, comment='#', fill_trailing_empty=True): 
	""" 
	return the iline-th line of the file which is non-empty and does not start with the comment ('#' by default).  
	if fill_trailing_emtpy is True then if iline is larger than the number of lines then return comma seperated empty values with the size the same as the header line. 

	iline could be either integer or slice instances, such as iline=slice(1, None, None) will return all lines after the first one. 

	Params
	------
	fn (str)
	iline=1 (int or slice instance)
	comment='#'
	fill_trailing_empty=True

	Return 
	------
	list of strings (lines)
	"""

	with open(fn, 'r') as f:
		data = f.read()
	lines = data.split('\n')

	lines_noncomment = []
	for line in lines:
		if len(line) > 0:
			if (line[0] != comment):
				lines_noncomment += [line]

	if iline < len(lines_noncomment) or isinstance(iline, slice): 
		return lines_noncomment[iline]

	elif fill_trailing_empty and len(lines_noncomment)>0:
		n_comma = lines_noncomment[0].count(',') 
		return "," * n_comma

	else:
		raise Exception("[batch] _extract_line_from_file iline exceeding the number of lines")

