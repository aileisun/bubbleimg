# filterzrange.py
# ALS 2016/05/02
"""
Determine the redshift range of a given line in a given filter. The outputs 
are .txt files: 

	filterboundary_0.2.txt
	filterboundary_0.6.txt
	filterboundary_0.8.txt
	HaNIIredshiftrange0.2.txt
	OIIIredshiftrange0.6.txt

"""
# from pylab import *

import filtertools
reload(filtertools)

import getzrange_batch
reload(getzrange_batch)

import getzrange_line
reload(getzrange_line)

def main():
	"""
	make the following files: 
	"""
	for threshold in [0.2, 0.6, 0.8]:
		filtertools.findFilterBounday(threshold=threshold,toplot=False)

	# write files
	getzrange_line.findzrange_wline_OIIIs(threshold=0.6)
	getzrange_line.findzrange_nline_HaNII(threshold=0.2)
	getzrange_line.findzrange_nline_HaNIISII(threshold=0.2)
	
	getzrange_batch.write_zranges()
		



if __name__ == '__main__':
	main()

