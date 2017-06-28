# test_filtertools.py

# WARNING this test set is not complete yet, not all functions in filtertools are covered

from .. import filtertools



def test_filtertools_getFilterBoundaries():

	band = 'r'
	survey = 'hsc'


	w1, w2 = filtertools.getFilterBoundaries(threshold=0.6, band=band, survey=survey, withunit=False)

	assert w1 == 5433
	assert w2 == 6967

	w1, w2 = filtertools.getFilterBoundaries(threshold=0.01, band=band, survey=survey, withunit=False)

	assert w1 == 5324
	assert w2 == 7070
