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

from . import filtertools
import imp
imp.reload(filtertools)

from . import getzrange_batch
imp.reload(getzrange_batch)

from . import getzrange_line
imp.reload(getzrange_line)

def main():
	"""
	make the following files: 
	"""

	for survey in ['sdss', 'hsc']: #, 'cfht', 'ukirt']:
		filtertools.writeFilterCentroids(survey=survey)
		filtertools.write_int_response_dlnl(survey=survey)
		filtertools.writeNormTransFunc(survey=survey)

		for threshold in [0.01, 0.2, 0.4, 0.6, 0.8]:
			filtertools.writeFilterBoundaries(threshold=threshold, toplot=False, survey=survey)

		# write files
		getzrange_line.findzrange_wline_OIIIs(threshold=0.6, survey=survey)
		getzrange_line.findzrange_wline_HaNII(threshold=0.6, survey=survey)
		getzrange_line.findzrange_nline_HaNII(threshold=0.2, survey=survey)
		getzrange_line.findzrange_nline_HaNII(threshold=0.4, survey=survey)
		getzrange_line.findzrange_nline_HaNIISII(threshold=0.2, survey=survey)
		getzrange_line.findzrange_nline_HaNIISII(threshold=0.4, survey=survey)
		getzrange_line.findzrange_nline_OII(threshold=0.2, survey=survey)
		getzrange_line.findzrange_nline_OII(threshold=0.4, survey=survey)
		getzrange_line.findzrange_nline_OIINeIII(threshold=0.2, survey=survey)
		getzrange_line.findzrange_nline_OIINeIII(threshold=0.4, survey=survey)
		getzrange_line.findzrange_nline_OIII(threshold=0.2, survey=survey)
		getzrange_line.findzrange_nline_OIII(threshold=0.4, survey=survey)

	
		getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='HaNIISII', wthreshold=0.6, nthreshold=0.2)
		getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='HaNII', wthreshold=0.6, nthreshold=0.4)

		filtertools.plotFilters2File(survey=survey)
		
	survey = 'hsc'
	# getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='OII', wthreshold=0.6, nthreshold=0.2)
	getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='OIINeIII', wthreshold=0.6, nthreshold=0.4)

	print(" wOIII nOII combination (usually i-r) has only been implemented on hsc but no other surveys")


if __name__ == '__main__':
	main()

