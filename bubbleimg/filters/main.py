# main.py
# ALS 2016/05/02
"""

Calculate the required files for spector from the filter transmission functions. 


Determine the redshift range of a given line in a given filter. The outputs 
are .txt files: 

sdss/ 
	filterboundary_0.01.txt
	filterboundary_0.2.txt
	filterboundary_0.4.txt
	filterboundary_0.6.txt
	filterboundary_0.8.txt
	filtercentroid.txt
	intresdlnl.txt
	zrange_nline_HaNII_0.2.txt
	zrange_nline_HaNII_0.4.txt
	zrange_nline_HaNIISII_0.2.txt
	zrange_nline_HaNIISII_0.4.txt
	zrange_nline_OII_0.2.txt
	zrange_nline_OII_0.4.txt
	zrange_nline_OIII_0.2.txt
	zrange_nline_OIII_0.4.txt
	zrange_nline_OIINeIII_0.2.txt
	zrange_nline_OIINeIII_0.4.txt
	zrange_wline_HaNII_0.6.txt
	zrange_wline_OIII_0.6.txt
	zranges_band_wOIII0.6_nHaNII0.4.txt
	zranges_band_wOIII0.6_nHaNIISII0.2.txt
	filters.pdf


hsc/

	filterboundary_0.01.txt
	filterboundary_0.2.txt
	filterboundary_0.4.txt
	filterboundary_0.6.txt
	filterboundary_0.8.txt
	filtercentroid.txt
	intresdlnl.txt
	zrange_nline_HaNII_0.2.txt
	zrange_nline_HaNII_0.4.txt
	zrange_nline_HaNIISII_0.2.txt
	zrange_nline_HaNIISII_0.4.txt
	zrange_nline_OII_0.2.txt
	zrange_nline_OII_0.4.txt
	zrange_nline_OIII_0.2.txt
	zrange_nline_OIII_0.4.txt
	zrange_nline_OIINeIII_0.2.txt
	zrange_nline_OIINeIII_0.4.txt
	zrange_wline_HaNII_0.6.txt
	zrange_wline_OIII_0.6.txt
	zranges_band_wOIII0.6_nHaNII0.4.txt
	zranges_band_wOIII0.6_nHaNIISII0.2.txt
	zranges_band_wOIII0.6_nOII0.2.txt
	zranges_band_wOIII0.6_nOIINeIII0.4.txt
	filters.pdf

wOIII nOII combination (usually i-r) has only been implemented on hsc but no other surveys

"""
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

	for survey in ['sdss', 'hsc']: 
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
	getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='OII', wthreshold=0.6, nthreshold=0.2)
	getzrange_batch.write_zranges(survey=survey, wline='OIII', nline='OIINeIII', wthreshold=0.6, nthreshold=0.4)



if __name__ == '__main__':
	main()

