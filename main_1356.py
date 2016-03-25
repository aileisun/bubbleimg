# main_1356.py
# ALS 2015/08/28

from astropy.table import Table
import os
# import alignstamp
# reload(alignstamp)

# import subtractimg
# reload(subtractimg)

import class_obsobj

import run_batch

def main_1356():

	"""
	PURPOSE: make image subtraction of J1356
	"""

	# setting
	dirout='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/test/SDSSJ1356+1026/'
	if not os.path.isdir(dirout): os.mkdir(dirout)

	# identifying object
	thelist=Table()
	thelist['RA']=[209.1921098]
	thelist['DEC']=[10.43585876]
	obj=class_obsobj.obsobj(thelist,catalog=None)

	# make img
	obj.dir_obj=dirout
	run_batch.objw_makeallimg(obj)

if __name__ == '__main__':
	main_1356()

