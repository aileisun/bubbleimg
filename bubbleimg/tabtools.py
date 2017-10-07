import io
import os
import numpy as np
import astropy.table as at

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
	with io.BytesIO() as f_temp: 
		row.write(f_temp, format='ascii.csv')
		rowstring = f_temp.getvalue()

	if not withheader:
		# take out header
		rowstring = '\n'.join(rowstring.splitlines()[1:]) + '\n'

	with open(fn, 'a') as f_to:
		f_to.write(rowstring)


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
		# print tab
		# print condi
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
	select_arr = [[str(tab[key][i]) == str(condi[key]) for i in range(len(tab))] for key in condi]
	select = np.all(select_arr, axis=0)

	# select_condis = np.array([tab[key] == condi[key] for key in condi])
	print select_arr
	# select = np.all(select_condis, axis=0)
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
	select = np.all(np.array([tab[key] == condi[key] for key in condi]), axis=0)
	if np.sum(select)>0:
		tab.remove_rows(select)

	return tab

