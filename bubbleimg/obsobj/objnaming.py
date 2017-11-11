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
		options: 'sdss', 'sdss_precise', 'j', 'j_precise'

	Return
	------
	name (str)
	"""
	if (ra is None) or (dec is None):
		return 'newobject'

	elif obj_naming_sys == 'sdss':
		J_str = getJstring_fromRaDec(ra, dec, precision='m')
		return 'SDSS'+J_str

	elif obj_naming_sys == 'sdss_precise':
		J_str = getJstring_fromRaDec(ra, dec, precision='s')
		return 'SDSS'+J_str

	elif obj_naming_sys == 'j':
		J_str = getJstring_fromRaDec(ra, dec, precision='m')
		return J_str

	elif obj_naming_sys == 'j_precise':
		J_str = getJstring_fromRaDec(ra, dec, precision='s')
		return J_str

	else: 
		raise ValueError("[objnaming] obj_naming_sys not recognized")


def getJstring_fromRaDec(ra, dec, precision='m'):
	""" 
	return J coordinate string given ra, dec 

	for example 'J1000+1242' 		if precision == 'm'
				'J100013+124226' 	if precision == 's'

	Params
	------
	ra (float)
	dec (float)
	precision='m' (str)
		'm' for minutes 's' for seconds
	"""
	c = SkyCoord(ra, dec, frame='icrs', unit='deg')

	if precision == 'm':
		ra_str = "%02d"%c.ra.hms[0] + "%02d"%c.ra.hms[1]
		dec_str = "%+03d"%c.dec.dms[0] + "%02d"%abs(c.dec.dms[1])
	elif precision == 's':
		ra_str = "%02d"%c.ra.hms[0] + "%02d"%c.ra.hms[1] + "%02d"%c.ra.hms[2]
		dec_str = "%+03d"%c.dec.dms[0] + "%02d"%abs(c.dec.dms[1]) + "%02d"%abs(c.dec.dms[2])
	else: 
		raise Exception('[objnaming] precision not understood')

	return 'J'+ra_str+dec_str

