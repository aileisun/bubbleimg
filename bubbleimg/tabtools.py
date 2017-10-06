import io
import os
import numpy as np
import astropy.table as at

def write_line(fn, line, condi, overwrite=False):
	"""
	write line (append) to file. If the line already exists, according to the condi conditions, then this line is overwritten (or not) depending on the overwrite parameter. 

	Params
	------
	fn (str)
	line (astropy tab)
	condi (dictionary)
		e.g., condi = {'imgtag': 'OIII5008_I'}
	overwrite=False
	"""
	if (not fn_has_line(fn, condi)) or overwrite:

		with io.BytesIO() as f_temp: 
			line.write(f_temp, format='ascii.csv')
			linestring = f_temp.getvalue()

		if os.path.isfile(fn):
			tab = at.Table.read(fn)
			tab = tab_delete_line(tab, condi)
			tab.write(fn, overwrite=True)
			# take out duplicated header
			linestring = '\n'.join(linestring.splitlines()[1:])

		with open(fn, 'a') as f_to:
			f_to.write(linestring)

	else: 
		print("[tabtools] skip writing line as it exists")


def fn_has_line(fn, condi):
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
		result = tab_has_line(tab, condi)

	else: 
		result = False

	return result


def tab_has_line(tab, condi):
	""" 
	return if table has a line with column (key) equals to value. 

	Params
	------
	tab: table
	condi (dictionary)
		{key: value, ...}, where key is the name of the column and value is the value that the column should take. 

		e.g., condi = {'imgtag': 'OIII5008_I'}
	"""
	select = np.all(np.array([tab[key] == condi[key] for key in condi]), axis=0)
	return np.sum(select) > 0



def fn_delete_line(fn, condi):
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
		tab = tab_delete_line(tab, condi)
		tab.write(fn, overwrite=True)


def tab_delete_line(tab, condi):
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

