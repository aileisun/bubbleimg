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
