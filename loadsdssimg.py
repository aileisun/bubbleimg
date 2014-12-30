# colordisp/loadsdssimg.py
# ipython -pylab

# purpose, load sdss image given list

from astropy.table import Table
from astroquery.sdss import SDSS




# testing ....

def loadsdssimg(list):
	crossid= Table.read('crossid_test.fits',format='fits')

	xid=Table(crossid[0])
	sp=SDSS.get_spectra(matches=Table(xid))
	im = SDSS.get_images(matches=xid, band='g')


	im[0].writeto('img_test.fits')


from astropy.coordinates import SkyCoord

c = SkyCoord(xid['ra'], xid['dec'], 'icrs', unit='deg')
result = SDSS.query_region(c,spectro=True,fields=['ra','dec','objid','run','rerun','camcol','field','z','plate','mjd','fiberID','specobjid','run2d','instrument','sciencePrimary'])

