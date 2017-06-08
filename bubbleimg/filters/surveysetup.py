# surveysetup.py
# ALS 2017/05/03

# """
# survey info
# """

import astropy.units as u

surveybands = {
                'sdss':['u','g','r','i','z'], 
                'hsc': ['g','r','i','z','y'],
                'ukirt': ['j', 'h', 'k'],
                'cfht': ['u']
                }

pixsize = {'sdss': 0.396 * u.arcsec,
           'hsc': 0.168 * u.arcsec
           }

waverange = {
			'sdss': [2980.*u.AA, 11230.*u.AA],
      'hsc': [3600.*u.AA, 11000.*u.AA],
			'cfht': [3000.*u.AA, 5540.*u.AA],
			}