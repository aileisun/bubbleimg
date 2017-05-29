# objnaming.py
# ALS 2017/05/26
"""
name object
"""

from astropy.coordinates import SkyCoord


def get_obj_name(ra, dec, obj_naming_sys='sdss'):
	"""
	return the name of the object given ra, dec based on the naming system

	Params
	------
	ra (float): in deg
	dec (float): in deg
	obj_naming_sys='sdss'

	Return
	------
	name (str)
	"""
	if obj_naming_sys == 'sdss':
		return getSDSSName_fromRADEC(ra, dec)
	else: 
		raise ValueError("[objnaming] obj_naming_sys not recognized")


def getSDSSName_fromRADEC(ra, dec):
	"""
	return SDSS name given 
	Params
	------
	ra (float) in deg
	dec (float) in deg
	"""
	c = SkyCoord(ra, dec, 'icrs', unit='deg')
	return getSDSSName_fromc(c)


def getSDSSName_fromc(c):
	return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+"%+03d"%c.dec.dms[0]+"%02d"%abs(c.dec.dms[1])
