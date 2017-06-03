# fromspec.py
# ALS 2015/08/10
"""
PURPOSE: Download and store SDSS spectrum of objects, 
		 calculate the continuum flux level at different bands, 
		 and return the ratio of continuum flux level between band1 and band2. 

USAGE:   getObjBandContiFLuxRatio(obj, band1='r', band2='z')

!!! WARNING !!! to be refactored into object oriented code

"""

import os
import numpy as np

from astropy.table import Table
from astroquery.sdss import SDSS
from astropy.io import fits
import astropy.units as u
from astropy import constants as const

import getcontspec
# reload(getcontspec)

from .. import filters

# import sys
# sys.path.append('../')
# import filters
# path_filters = os.path.dirname(filters.__file__)+'/'

def obj_measure_contiABmags(obj, bands=['u','g','r','i','z']):
	"""
	PURPOSE: measure the AB magnitudes the spec continuum in each band and return as a table
	"""

	tabout=Table()

	for band in bands:
		tabout['contiMag_'+band]=[getObjBandContiABmag(obj, band=band)]
	return tabout


def getObjBandContiABmag(obj, band='r', survey='sdss'):
	"""
	PURPOSE: get the AB magnitude of the SDSS spec continuum in a band

	PAREMATERS: 
		obj (obsobj) : obsobj instance
		band='r' (string)

	RETURN m  AB magnitude

	DESCRIPTION: 
		m_AB = -2.5 * log_10(f_nu/ (erg/s/cm2/Hz) ) - 48.600
		f_nu = lambda^2 / c * f_lambda
	"""	
	# set up
	u_fnu = u.Unit('erg/cm^2/s/Hz')

	fl = getObjBandContiFluxDensity_dEdl(obj, band=band, survey=survey, wunit=True)

	wave_ctr = filters.filtertools.getFilterCentroids(band=band, survey=survey, withunit=True)
	fnu = (wave_ctr**2/const.c*fl).to(u_fnu)

	m = -2.5 * np.log10(fnu/u_fnu) - 48.600

	return m


def getObjBandContiRatio_dEdnu(obj, band1='r', band2='z', survey='sdss'):
	"""
	PURPOSE: get the ratio of continuum flux density (dE/dnu) between two bands. 

	DESCRIPTION: 
		Using getObjBandContiRatio_dEdl() to get ratio and transform it from dE/dl to dE/dnu

			dE/dnu = - (lambda^2/c) * dE/dlambda

			ratio12_dEdnu = (lambda1/lambda2)^2 * ratio12_dEdl

		For details, see getObjBandContiRatio_dEdl()
		This ratio can be applied to SDSS images in units of 'nanomaggies' (dE/dnu)
	"""

	ratio12_dEdl = getObjBandContiRatio_dEdl(obj, band1=band1, band2=band2, survey=survey)

	w1 = filters.getFilterCentroids(survey=survey, band=band1, withunit=True)
	w2 = filters.getFilterCentroids(survey=survey, band=band2, withunit=True)

	ratio12_dEdnu = (w1/w2)**2 * ratio12_dEdl
	return ratio12_dEdnu.value


def getObjBandContiRatio_dEdl(obj, band1='r', band2='z', survey='sdss'):
	"""
	PURPOSE: get the ratio of continuum flux density (dE/dlambda) between two bands. 

	PAREMATERS: 
		obj (obsobj) : obsobj instance
		band1='r' (string)
		band2='z' (string)

	RETURN f1/f2 (float)

	DESCRIPTION: 
		using function getcontspec.getMedianFilteredConti() to get continuum levels
					   convolveSpecWFilter() to convolve with filter
		as the spectrum is in units of '1E-17 erg/cm^2/s/Ang', the ratio is of dE/dlambda
	"""
	# get flux levels
	fl1=getObjBandContiFluxDensity_dEdl(obj, band=band1, survey=survey, wunit=False)
	fl2=getObjBandContiFluxDensity_dEdl(obj, band=band2, survey=survey, wunit=False)

	# return ratio f1/f2
	return fl1/fl2


def getObjBandContiFluxDensity_dEdl(obj, band='r', survey='sdss', wunit=False):
	"""
	PURPOSE: get the SDSS spec continuum flux density (dE/dlambda) in a band

	PAREMATERS: 
		obj (obsobj) : obsobj instance
		band='r' (string)

	RETURN f (flux density quantity in units of '1E-17 erg/cm^2/s/Ang')

	DESCRIPTION: 
		using function getcontspec.getMedianFilteredConti() to get continuum levels
					   convolveSpecWFilter() to convolve with filter
		flux density in units of '1E-17 erg/cm^2/s/Ang'
	"""

	# load SDSS spectrum
	spec, lcoord=obj.sdss.get_speclcoord(wunit=wunit)
	try: lcoord=lcoord.value
	except: pass

	# get continuum
	speccont, lcoordcont = getcontspec.getMedianFilteredConti(spec, lcoord, obj.sdss.z, toplot=False)

	# get filter function
	specfilter, lcoordfilter = filters.getFilterResponseFunc(band=band, survey=survey)

	# convolve to get flux density
	fl = convolveSpecWFilter(speccont, lcoordcont ,specfilter, lcoordfilter)

	# return flux density
	return fl


def getObjBandABmag(obj, band='r', survey='sdss'):
	"""
	PURPOSE: get the AB magnitude of the SDSS spec in a band

	PAREMATERS: 
		obj (obsobj) : obsobj instance
		band='r' (string)

	RETURN m  AB magnitude

	DESCRIPTION: 
		m_AB = -2.5 * log_10(f_nu/ (erg/s/cm2/Hz) ) - 48.600
		f_nu = lambda^2 / c * f_lambda
	"""	
	# set up
	u_fnu = u.Unit('erg/cm^2/s/Hz')

	fl = getObjBandFluxDensity_dEdl(obj, band=band, survey=survey, wunit=True)

	wave_ctr = filters.filtertools.getFilterCentroids(band=band, survey=survey, withunit=True)
	fnu = (wave_ctr**2/const.c*fl).to(u_fnu)

	m = -2.5 * np.log10(fnu/u_fnu) - 48.600

	return m


def getObjBandFluxDensity_dEdl(obj, band='r', survey='sdss', wunit=False):
	"""
	PURPOSE: get the SDSS spec flux density (dE/dlambda) in a band

	PAREMATERS: 
		obj (obsobj) : obsobj instance
		band='r' (string)

	RETURN f (flux density quantity in units of '1E-17 erg/cm^2/s/Ang')

	DESCRIPTION: 
		using function convolveSpecWFilter() to convolve with filter
		flux density in units of '1E-17 erg/cm^2/s/Ang'
	"""

	# load SDSS spectrum
	spec, lcoord=obj.sdss.get_speclcoord(wunit=wunit)
	try: lcoord=lcoord.value
	except: pass

	# get filter function
	specfilter, lcoordfilter = filters.getFilterResponseFunc(band=band, survey=survey)

	# convolve to get flux density
	fl = convolveSpecWFilter(spec, lcoord ,specfilter, lcoordfilter)

	# return flux density
	return fl


def R_lambda(wavelength, band='r', survey='sdss'):
	"""
	PURPOSE: Interpolated filter response function as a function of wavelength [AA] given band
	PAREMATERS: 
		wavelength (float or array) [AA]
		band='r'   (string)

	RETURN: reponse funciton (not normalized) at the wavelength(s)
	"""
	from scipy.interpolate import interp1d
	# read filter response fuction
	specfilter, lcoordfilter=filters.getFilterResponseFunc(band=band, survey=survey)
	f = interp1d(lcoordfilter, specfilter,kind='linear',bounds_error=False,fill_value=0.)
	return f(wavelength)



def convolveSpecWFilter(spec, lcoord ,specfilter, lcoordfilter):
	"""
	PUPROSE: convolve spectrum with (normalized) filter response funciton
			 and return a filter response funciton weighted averaged spectrum intenisty

	DESCRIPTION: reponse funciton is linearly interpolated to each spec data points to get the weight. 

	RETURN: specavg (float)
	"""
	from scipy.interpolate import interp1d
	# spline response function
	f = interp1d(lcoordfilter, specfilter,kind='linear',bounds_error=False,fill_value=0.)

	specavg=np.average(spec,weights=f(lcoord))

	# glue back the unit if spec has one
	try: spec.unit
	except: pass
	else: specavg=specavg*spec.unit		

	return specavg



def plotSpecwFilters(obj, survey='sdss'):
	"""

	PURPOSE: Overplot filter response function on top of spectrum, to see what lines are in what bands. 

	PAREMATERS:    obj

	WRITE OUTPUTS: obj.dir_obj+'spec.pdf'
	"""
	import matplotlib.pyplot as plt
	# settings
	fileout=obj.dir_obj+'spec.pdf'

	# load SDSS spectrum
	spec, lcoord=obj.sdss.get_speclcoord()

	# plotting
	plt.figure(1,figsize=(12, 6))
	plt.clf()
	plt.plot(0.,0.,ls='',color='black',label='z='+'%.3f'%obj.sdss.z)
	plt.plot(lcoord,spec/max(spec),color='black',lw=2)

	for band in filters.filtertools.surveybands[survey]:
		specfilter, lcoordfilter=filters.getFilterResponseFunc(band=band, survey=survey)
		plt.plot(lcoordfilter, specfilter/max(specfilter), label=band)

	plt.legend(loc='upper right')
	plt.xlim(2980.0, 11230.0)
	plt.ylim(0., 1.)
	plt.savefig(fileout)


def getSDSSspeclineratio(obj,line1='H_beta',line2='[O_III] 5007'):
	"""
	PURPOSE: Get SDSS fit line ratio from SDSS file spec-PLATE-MJD-FIBER.fits HDU[3]. 
	PARAMETERS:  obj,
				 line1='H_beta',
				 line2='[O_III] 5007' (other line options: e.g. '[O_III] 4959')
	RETURN: 	 r (float) line ratio

	DESCRIPTION: the flux used is the area under the best fit Gaussian, so may differ from true flux if line shape deviates from Gaussian. 
	"""
	table=obj.sdss.get_spec()[3].data
	f1=table[table['LINENAME']==line1]['LINEAREA'][0]
	f2=table[table['LINENAME']==line2]['LINEAREA'][0]

	return f1/f2


