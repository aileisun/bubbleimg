# spector.py 
# ALS 2017/06/01

import sys
import numpy as np
import astropy.table as at
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
import astropy.constants as const
import modelBC03
import os

from ..obsobj import Operator
from .. import filters
from . import getconti
from . import extrap
from . import linelist
from . import lineflux


class Spector(Operator):

	def __init__(self, **kwargs):
		"""
		Spector, an operator on spectrum

		Params
		----------
		Operator params:
			/either
				obj (object of class obsobj): with attributes ra, dec, dir_obj
			/or  
				ra (float)
				dec (float)
				/either
					dir_obj (string)
				/or 
					dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'

		survey (str): 
			survey of the photometric system
			if not provided, use self.obj.survey. Raise exception if self.obj.survey does not exist. 
		z (float):  
			redshift, if not provided, use obj.z or obj.sdss.z
		survey_spec (str): 
			survey of the spectrum
			if not provided, use the instrument of sdss_xid.csv
		decompose_method = 'modelBC03' (str)
			'modelBC03' or 'running_median'
			the method for decomposing spectrum into continuum and lines, and extrapolate the continum. 


		Attributes
		----------
		Operator Attributes:	
			obj (instance of objObj)
			ra (float)
			dec (float)
			dir_obj (string)

		survey (str): e.g., 'hsc'
			survey of the photometric system
		survey_spec (str): e.g., 'sdss' or 'boss'
			survey of the spectrum
		bands (list): 
			e.g., ['g', 'r', 'i', 'z', 'y'] for survey = 'hsc'
		z (float):  
		decompose_method (str)
			'modelBC03' or 'running_median'
		conti_model
			if decompose_method is modelBC03 then this is an instance of modelBC03
		"""
		
		super(Spector, self).__init__(**kwargs)

		# set survey
		if hasattr(self.obj, 'survey'):
			default_survey = self.obj.survey
			self.survey = kwargs.pop('survey', default_survey)
		else: 
			self.survey = kwargs.pop('survey')

		# set survey_spec
		self.survey_spec = kwargs.pop('survey_spec', 'auto')

		if self.survey_spec in ['sdss', 'boss', 'eboss', 'auto']:
			if self.survey_spec == 'auto':
				self.obj.add_sdss(toload_photoobj=False)
				self.survey_spec = self.obj.sdss.instrument.lower()

			# # sanity check - spec name consistent
			# if (self.survey_spec != self.obj.sdss.instrument.lower()):
			# 	raise Exception("[spector] survey_spec inconsistent with sdss_xid.csv:instrument")

		# set z
		if 'z' in kwargs:
			self.z = kwargs.pop('z')

		elif hasattr(self.obj, 'z'):
			self.z = self.obj.z

		elif self.survey_spec in ['sdss', 'boss', 'eboss', 'auto']:
			self.obj.add_sdss(toload_photoobj=False)
			self.z = kwargs.pop('z', self.obj.sdss.z) 

		# set self.decompose_method
		self.decompose_method = kwargs.pop('decompose_method', 'modelBC03')
		self.conti_model = None

		# set others
		self.bands = filters.filtertools.surveybands[self.survey]
		self.waverange = filters.filtertools.waverange[self.survey]

		# define paths
		self.fp_spec = self.dir_obj+'spec.fits'
		self.fp_spec_decomposed = self.dir_obj+'spec_decomposed.ecsv'
		self.fp_spec_contextrp = self.dir_obj+'spec_contextrp.ecsv'
		self.fp_spec_mag = self.dir_obj+'spec_mag.csv'
		self.fp_spec_lineflux = self.dir_obj+'spec_lineflux.csv'
		self.fp_spec_linefrac = self.dir_obj+'spec_linefrac.csv'

		self.spec, self.ws = self.get_spec_ws()
		self.u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)') # default unit of spec
		self.u_ws = u.AA # default unit of ws


	def get_spec_ws(self, forceload_from_fits=False):
		"""
		read spec and ws, either from spec_decomposed.ecsv or spec.fits (if forced or spec_decomposed.ecsv does not exist)

		Param
		-----
		forceload_from_fits=False:
			if true, then load from fits. 

		Return
		------
		spec (astropy table column with unit)
		ws (astropy table column with unit)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_ws = u.AA
		"""

		fn = self.fp_spec_decomposed

		if not os.path.isfile(fn) or forceload_from_fits:
			spec, ws, __ = self.__read_spec_ws_ivar_from_fits()
			return spec, ws
		else: 
			tab = at.Table.read(fn, format='ascii.ecsv')
			return tab['spec'], tab['ws']


	def get_spec_ws_from_spectab(self, component, fn=None):
		"""
		read certain spec table.ecsv to get spec of a given component
		if fn not specified then use the default file that would contain the specified component,
		e.g.:
			'spec_decomposed.ecsv' for component = "all", "cont", "line"
			'spec_contextrp.ecsv'    for component = "contextrp"

		Param
		-----
		component: 
			either ['all', 'line', 'cont', 'contextrp']
		fn='spec_decomposed.ecsv'

		Return
		------
		spec (astropy table column with unit)
		ws (astropy table column with unit)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_ws = u.AA
		"""
		if fn == None:
			if component in ['all', 'cont', 'line']:
				fn = self.fp_spec_decomposed
				self.make_spec_decomposed_ecsv(overwrite=False)
			elif component in ['contextrp']:
				fn = self.fp_spec_contextrp
				self.make_spec_contextrp_ecsv(overwrite=False)
			else: 
				raise Exception("[spector] component not recognized")

		tab = at.Table.read(fn, format='ascii.ecsv')
		col = self.__get_spectab_colname(component)

		return tab[col], tab['ws']


	def make_spec_decomposed_ecsv(self, overwrite=False):
		""" 
		saving seperated continuum and line spectrum in csv file 

		Params
		------
		self
		overwrite=False

		Return
		------
		status

		note: this part of the code can be refactorized better. 

		"""

		fn = self.fp_spec_decomposed

		if (not os.path.isfile(fn)) or overwrite:

			spec, ws = self.get_spec_ws()

			ws_uless = np.array((ws/self.u_ws).to(u.dimensionless_unscaled))
			spec_uless = np.array((spec/self.u_spec).to(u.dimensionless_unscaled))

			iscon, speccont, specline, __, model = getconti.decompose_cont_line_t2AGN(spec_uless, ws_uless, self.z, method=self.decompose_method)

			self.conti_model = model # modelBC03 instance if set method='modelBC03', otherwise None.

			tab = at.Table([ws_uless, spec_uless, speccont, specline, iscon], names=['ws', 'spec', 'speccont', 'specline', 'iscon'])
			tab['ws'].unit = self.u_ws
			tab['spec'].unit = self.u_spec
			tab['speccont'].unit = self.u_spec
			tab['specline'].unit = self.u_spec

			tab.write(fn, format='ascii.ecsv', overwrite=overwrite)

			# sanity check: units are identical
			units = [tab[col].unit for col in ['spec', 'speccont', 'specline']]
			if len(set(units)) > 1:
				raise Exception("[spector] units in table spec_decomposed are not identical")

		status = os.path.isfile(fn)
		return status 


	def make_spec_contextrp_ecsv(self, overwrite=False, refit=False):
		""" 
		extrapolate continuum to cover all of the wavelength range of filters 
		there are two methods:
			for self.conti_model modelBC03: use the bestfit
			for running_median: polynomial fit

		"""
		fn = self.fp_spec_contextrp

		if (not os.path.isfile(fn)) or overwrite:
			speccont, ws = self.get_spec_ws_from_spectab(component='cont')
			ws_uless = np.array((ws/self.u_ws).to(u.dimensionless_unscaled))
			speccont_uless = np.array((speccont/self.u_spec).to(u.dimensionless_unscaled))

			l0, l1 = self.waverange

			if self.decompose_method == 'modelBC03':
				if self.conti_model is None or refit:
					m = modelBC03.modelBC03(extinction_law='none')
					m.fit(ws=ws_uless, spec=speccont_uless, z=self.z)
				else: 
					m = self.conti_model # reuse 
				speccon_ext = m.bestfit
				ws_ext = m.ws_bestfit

			elif self.decompose_method == 'running_median':			
				speccon_ext, ws_ext = extrap.extrap_to_ends(ys=speccont_uless, xs=ws_uless, x_end0=l0, x_end1=l1, polydeg=1, extbase_length=2000.)
			else:
				raise Exception("method is not recognized")

			assert (np.min(ws_ext*self.u_ws) < l0) & (np.max(ws_ext*self.u_ws) > l1)
			col_contextrp = self.__get_spectab_colname('contextrp')
			tab = at.Table([ws_ext, speccon_ext], names=['ws', col_contextrp])
			tab['ws'].unit = self.u_ws
			tab[col_contextrp].unit = self.u_spec
			tab.write(fn, format='ascii.ecsv', overwrite=overwrite)

		status = os.path.isfile(fn)
		return status 


	def make_spec_mag(self, overwrite=False):
		"""
		make table spec_mag.csv that contains the convolved spectral magnitude and fnu in each band

		Params
		------
		self
		overwrite=False

		Return
		------
		status
		"""
		#==========================================================================

		fn = self.fp_spec_mag

		self.make_spec_decomposed_ecsv(overwrite=False)

		if not os.path.isfile(fn) or overwrite:
			print("[spector] making spec_mag")
			tabmag = at.Table()
			tabfnu = at.Table()

			for component in ['all', 'cont', 'line', 'contextrp']:
				for band in self.bands:
					colfnu = self.__get_specmag_colname(band, component=component, fluxquantity='fnu')
					colmag = self.__get_specmag_colname(band, component=component, fluxquantity='mag')
					try:
						fnu = self._calc_Fnu_in_band(band=band, component=component)
					except KeyboardInterrupt:
						sys.exit(0) 
					except:
						print(("[spector] skip calculating fnu of {} in band {}".format(component, band)))
					else: 
						mag = fnu.to(u.ABmag)
						fnu_nm = fnu.to(u.nanomaggy)				
						tabmag[colmag] = [mag.value]
						tabfnu[colfnu] = [fnu_nm.value]

			tab = at.hstack([tabmag, tabfnu])
			tab.meta['comments'] = [
									"survey_photo: {}".format(self.survey),
									"survey_spec: {}".format(self.survey_spec),
									"unit_mag: ABmag",
									"unit_fnu: nanomaggy",
									]

			tab.write(fn, comment='#', format='ascii.csv', overwrite=overwrite)
		else:
			print("[spector] skip making spec_mag as file exists")

		status = os.path.isfile(fn)
		return status 



	def get_spec_mag_tab(self):
		""" return spec_mag table"""
		self.make_spec_mag(overwrite=False)
		return at.Table.read(self.fp_spec_mag, comment='#', format='ascii.csv')


	def get_spec_mag_value(self, band, component='all', fluxquantity='mag'):
		""" 
		return requested value, reading from spec_mag.csv

		Params
		------
		component='line':
			or 'all', 'cont', 'contextrp'
		fluxquantity='mag'
		band='i'

		Return
		------
		x (float): the value
		"""
		tab = self.get_spec_mag_tab()
		col = self.__get_specmag_colname(band=band, component=component, fluxquantity=fluxquantity)
		
		return tab[col][0]


	def get_fnu_ratio_band1_over_band2(self, band1, band2, component='all'):
		""" 
		return fnu_band1/fnu_band2

		Params
		------
		band1, band2, component='all'

		Return
		------
		x (float): the ratio
		"""
		fluxquantity = 'fnu'
		tab = self.get_spec_mag_tab()
		col1 = self.__get_specmag_colname(band=band1, component=component, fluxquantity=fluxquantity)
		col2 = self.__get_specmag_colname(band=band2, component=component, fluxquantity=fluxquantity)

		ratio = tab[col1][0]/tab[col2][0]
		
		return ratio


	def plot_spec(self, wfilters=True, wconti=True, wcontextrp=True, wline=True, wspec=True, overwrite=True):
		"""
		Params
		------
		self
		wfilters=True
		wconti=False

		Return
		------
		status (bool)
		"""
		fn = self.dir_obj+'spec.pdf'

		if not os.path.isfile(fn) or overwrite:
			plt.close('all')
			plt.figure(1, figsize=(12, 6))
			plt.clf()


			plt.plot(0., 0., ls='', color='black', label='z='+'%.3f'%self.z) # show z in legend

			if wfilters:
				for band in self.bands:
					trans, ws_trans = self._get_norm_trans_func(band=band)
					plt.plot(ws_trans, trans/max(trans), label=band)

			if wspec:
				spec, ws = self.get_spec_ws()
				norm = max(spec)			
				plt.plot(ws, spec/norm, color='0.3', lw=1.5, label='__nolabel__')

			if wline:
				speccont, ws = self.get_spec_ws_from_spectab(component='cont')
				specline, ws = self.get_spec_ws_from_spectab(component='line')
				specplot = specline+speccont
				specplot[specline==0] = np.nan
				plt.plot(ws, specplot/norm, color='black', lw=1.5, label='line')

			if wconti:
				speccont, ws = self.get_spec_ws_from_spectab(component='cont')
				plt.plot(ws, speccont/norm, color='cyan', lw=2, label='continuum')

			if wcontextrp:
				speccextrp, ws = self.get_spec_ws_from_spectab(component='contextrp')
				plt.plot(ws, speccextrp/norm, color='grey', lw=0.5, label='extrapolated conti')


			plt.legend(loc='upper right')
			if self.survey_spec == 'sdss':
				plt.xlim(2980.0, 11230.0)
			elif self.survey_spec == 'boss':
				plt.xlim(3500.0, 12000.0)
			plt.ylim(0., 1.)
			plt.savefig(fn, overwrite=overwrite)

		status = os.path.isfile(fn)
		return status 


	def calc_fline_over_fnuband(self, band, line):
		""" 
		calculate the conversion ratio r = (flux_line / fnu_band), to convert continuum subtracted image to line intensity map

		by multiplying this ratio to the band image (in fnu [erg/s/cm^2/Hz]) one can obtain the observed flux of the line (in f [erg/s/cm^2]). Notice that this is flux in observed frame, for measurements one still needs to transform it to the rest frame (depending on z). 

		r = (c / (T(w_k) * w_k)) * frac_k  := dnu * frac_k

			where 
				dnu    =  (c / (T(w_k) * w_k))
				frac_k = (f_k * T(w_k) * w_k) / sum(f_i * T(w_i) * w_i)

		Params
		------
		self
		band (str)
		line (str)

		Return
		------
		ratio (quantity in unit of Hz)
		"""
		f, __ = self._get_line_flux(line=line, wunit=False)
		w = self._get_line_obs_wave(line=line, wunit=False)
		T = self._get_norm_trans(wavelength=w, band=band, bounds_error=True)

		dnu = const.c / (T * w * u.AA)

		frac = self._get_line_frac(band=band, line=line)

		# sanity check: all strong lines are considered
		tab_linefrac = at.Table.read(self.fp_spec_linefrac, format='ascii.csv', comment='#')
		stronglines = self._list_stronglines_in_band(band=band)
		for line in stronglines:
			if 'frac_{}'.format(line) not in tab_linefrac.colnames:
				raise Exception("[spector] strong line {} in band is not contained in spec_linefrac.csv".format(line))

		r = (dnu * frac).to(u.Hz)
		return r


	def _get_line_frac(self, band, line='OIII5008'):
		""" 
		read line flux from file spec_linefrac.csv. 

		Params
		------
		band, line='OIII5008', wunit=False

		Return
		------
		frac (float)
		"""
		fn = self.fp_spec_linefrac

		if not os.path.isfile(fn):
			self.make_linefrac(band=band, overwrite=False)

		tab = at.Table.read(fn, format='ascii.csv', comment='#')

		if tab['lineband'][0] != band:
			raise Exception("[spector] _get_line_frac the band required is not provided by the current spec_linefrac.csv")

		linetag = 'frac_{}'.format(line)
		frac = tab[linetag][0]

		return frac


	def make_linefrac(self, band, lines=['NeIII3870', 'NeIII3969', 'Hg', 'Hb', 'OIII4960', 'OIII5008', 'OI6302', 'OI6366'], tofixOIIIratio=True, overwrite=False):
		"""
		make file spec_linefrac.csv that contains the fraction each of the strong lines have in a specific band.
		Columns: f_{band}_{line}, T_{band}_{line}, w_{band}_{line}, frac_{band}_{line}

		The fraction is based on the f*T*w of the line. Only the strong lines are listed. 

		If tofixOIIIratio = True, then the ratio between OIII5008 and OIII4960 is fixed to the theoretical ratio of 2.98, see 
		Storey + 2000. http://adsabs.harvard.edu/abs/2000MNRAS.312..813S. 

		Params
		------
		band
		lines=['NeIII3870', 'NeIII3969', 'Hg', 'Hb', 'OIII4960', 'OIII5008', 'OI6302', 'OI6366']
		tofixOIIIratio=True
		overwrite=False

		Return
		------
		status
		"""
		fn = self.fp_spec_linefrac

		self.make_lineflux(overwrite=overwrite)

		if not os.path.isfile(fn) or overwrite:
			print("[spector] making spec_linefrac")

			tab = at.Table([[band]], names=['lineband'])

			fwt_sum = 0.
			for line in lines:
				f, __ = self._get_line_flux(line=line, wunit=False)
				w = self._get_line_obs_wave(line=line, wunit=False)
				T = self._get_norm_trans(wavelength=w, band=band, bounds_error=False)

				fwt = max(f*w*T, 0)
				fwt_sum = fwt_sum + fwt

				col_new = at.Table([[f], [w], [T], [fwt]], names=['f_{}'.format(line), 'w_{}'.format(line), 't_{}'.format(line), 'fwt_{}'.format(line)])
				tab = at.hstack([tab, col_new])

			for line in lines:
				frac = tab['fwt_{}'.format(line)][0] / fwt_sum
				col_new = at.Table([[frac]], names=['frac_{}'.format(line)])
				tab = at.hstack([tab, col_new])

			if tofixOIIIratio:
				r = 2.98
				frac_OIIItotal = tab['frac_OIII4960'] + tab['frac_OIII5008']
				tab['frac_OIII5008'] = frac_OIIItotal * r / (1.+r)
				tab['frac_OIII4960'] = frac_OIIItotal * 1. / (1.+r)

			tab.write(fn, format='ascii.csv', overwrite=overwrite)

		else:
			print("[spector] skip making spec_linefrac as file exists")

		status = os.path.isfile(fn)
		return status 


	def make_lineflux(self, lines=['NeIII3870', 'NeIII3969', 'Hg', 'Hb', 'OIII4960', 'OIII5008', 'OI6302', 'OI6366'], u_flux=u.Unit("1E-17 erg cm-2 s-1"), overwrite=False):
		""" 
		make file spec_lineflux.csv that contains the flux of the specified lines. The fluxes are calculated by integrating the line component of the spectrum over a window of +/- 1400 km/s. 

		WARNING: currently only Hb, and OIII lines are supported. For lines that are overlapped, e.g., Ha and NII, the current implemenation will double count the flux. 

		Params
		------
		lines=['Hb', 'OIII4960', 'OIII5008']
		u_flux=u.Unit("1E-17 erg cm-2 s-1")
			the unit of the output
		overwrite=False

		Return
		------
		status (bool)
		"""
		fn = self.fp_spec_lineflux

		self.make_spec_decomposed_ecsv(overwrite=False)

		if not os.path.isfile(fn) or overwrite:
			print("[spector] making spec_lineflux")
			tab = at.Table()

			for line in lines:
				f, ferr = self._calc_line_flux(line=line, u_flux=u_flux, wunit=False)

				col_new = at.Table([[f], [ferr]], names=['f_{}'.format(line), 'ferr_{}'.format(line)])
				tab = at.hstack([tab, col_new])

			tab.meta['comments'] = ["unit_flux: {}".format(u_flux.to_string()),]

			tab.write(fn, comment='#', format='ascii.csv', overwrite=overwrite)

		else:
			print("[spector] skip making spec_lineflux as file exists")

		status = os.path.isfile(fn)
		return status 


	def _calc_line_flux(self, line='OIII5008', dv=1400*u.km/u.s, u_flux=u.Unit("1E-17 erg cm-2 s-1"), wunit=False):
		"""
		calculate the flux of the line by trpz integrating the decomposed emission line component of the spectrum over a range of +/- dv 

		WARNING: currently only Hb, and OIII lines are supported. For lines that are overlapped, e.g., Ha and NII, the current implemenation will double count the flux. 

		WARNING: the errors of the flux is propogated from the variance ('ivar' column) of the sdss spectrum. However, sdss's noise could be correlated between neighbors but such a covariance is not quantified, which will result in 10-20% in the error estimates.  We here artifically boost the flux error by 20% to incoporate such uncertainty on the error, see:
		http://www.sdss.org/dr12/spectro/caveats/

		Params
		------
		line='OIII5008' (str)

		dv=1400*u.km/u.s (quantity)
			half of the width over which the spectrum is integrated to calculate flux

		u_flux=u.Unit("1E-17 erg cm-2 s-1")
			the unit of the output

		wunit=False
			whether the output to come with unit or not

		Return
		------
		f (float or quantity)
			unit is in the dimension of erg cm-2 s-1
		ferr (float or quantity)
			same unit as f
		"""

		# sanity check
		if line not in ['NeIII3870', 'NeIII3969', 'Hg', 'Hb', 'OIII4960', 'OIII5008', 'OI6302', 'OI6366']:
			raise Exception("[spector] _calc_line_flux does not support lines other than Hb and OIII as those are not tested. ")

		# get w range
		beta = (dv/const.c).to_value(u.dimensionless_unscaled)
		w = self._get_line_obs_wave(line=line, wunit=False)
		w0 = w*(1-beta)
		w1 = w*(1+beta)

		# get spectrum
		spec, ws = self.get_spec_ws_from_spectab(component='line')
		__, __, ivar = self.__read_spec_ws_ivar_from_fits()

		f, ferr = lineflux.calc_line_flux(spec, ws, ivar, w0, w1, u_flux)

		# artificially boost the error to account for pixel covariance
		ferr = ferr * 1.2

		if wunit:
			return f, ferr
		else: 
			return f.to_value(u_flux), ferr.to_value(u_flux)


	def _get_line_flux(self, line='OIII5008', wunit=False):
		""" 
		read line flux from file spec_lineflux.csv. For details, see make_spec_lineflux(). 

		Params
		------
		line='OIII5008', wunit=False

		Return
		------
		f, ferr
		"""
		fn = self.fp_spec_lineflux

		if not os.path.isfile(fn):
			self.make_lineflux(overwrite=False)

		tab = at.Table.read(fn, format='ascii.csv', comment='#')

		f = tab['f_{}'.format(line)][0]

		ferr = tab['ferr_{}'.format(line)][0]

		if not wunit:
			return f, ferr

		else:
			u_flux = u.Quantity(tab.meta['comments'][0].split(': ')[1])
			return f*u_flux, ferr*u_flux


	def _get_line_flux_sdss(self, line='OIII5008', wunit=False): 
		""" 
		Flux of line from SDSS spec.fits higher extensions 
		It is measured by Gaussian fit, for details, see:
		http://classic.sdss.org/dr7/dm/flatFiles/spZline.html

		PARAMS
		------
		line = 'OIII5008' (str)
		wunit = False

		Return
		------
		f (float): 
			flux
			or astropy quantify with units u.Unit("1E-17 erg cm-2 s-1")
		ferr (float)
			error on flux
			or astropy quantify with units u.Unit("1E-17 erg cm-2 s-1")
		"""

		# get sdss linename 
		linename = linelist.sdssLINENAME[line]

		# read flux
		table = self.obj.sdss.get_spec()[3].data
		i = [table['LINENAME']==linename]
		f = table[i]['LINEAREA'][0]
		ferr = table[i]['LINEAREA_ERR'][0]

		if not wunit:
			return f, ferr
		else:
			return f*u.Unit("1E-17 erg cm-2 s-1"), ferr*u.Unit("1E-17 erg cm-2 s-1")


	def _get_line_obs_wave(self, line='OIII5008', wunit=False): 
		"""
		PARAMS
		------
		line = 'OIII5008' (str)
		wunit = False

		Return
		------
		w (float):
			redshifted wavelength of line
			or astropy quantify with units u.Unit("AA")
		"""
		w = filters.getllambda(ion=line, vacuum=True) * (1. + self.z)

		try: len(w)
		except:
			pass
		else:
			if len(w) == 1:
				w = w[0]
			else:
				raise Exception('got more than one wavelength')

		if not wunit:
			return w
		else:
			return w*u.Unit("AA")


	def _list_stronglines_in_band(self, band, threshold=0.01):
		""" 
		return a list of strong lines in band where the filter transmission function is higher than threshold (in fraction) 

		Params
		------
		self
		band (str)
		threshold (float):
			the fractional threshold (relative to the filter peak) that defines the wavelength boundary of a band
		wunit (bool)

		Return
		------
		llist (array of str)
			e.g., ['OIII5008', 'OIII4960', 'Hb', ...]
		"""
		w1, w2 = filters.filtertools.getFilterBoundaries(threshold=threshold, band=band, survey=self.survey, withunit=False)

		llist = []
		for line in linelist.strongline:
			w = self._get_line_obs_wave(line=line, wunit=False)
			if (w > w1) & (w < w2):
				llist += [line]

		return llist


	def __read_spec_ws_ivar_from_fits(self, u_spec=1.e-17*u.Unit('erg / (Angstrom cm2 s)'), u_ws=u.AA, wunit=True):
		"""
		Params
		------
		self
		u_spec=1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_ws=u.AA
		wunit=True

		Return
		------
		spec (nparray)
		ws (nparray)
		ivar (nparray)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_ws = u.AA
		"""
		fn = self.fp_spec
		if self.survey_spec in ['sdss', 'boss', 'boss']:
			if os.path.isfile(fn):
				hdus = fits.open(fn)
			else: 
				raise IOError("[Spector] spec fits file does not exist")

			hdus[0].header
			spectable = hdus[1].data
			spec, ws, ivar = spectable['flux'], 10.**spectable['loglam'], spectable['ivar']

			if wunit:
				spec = spec*u_spec
				ws = ws*u_ws
				ivar = ivar / (u_spec**2)

			# instrument_header = at.Table(hdus[2].data)['INSTRUMENT'][0].lower()
			# if self.survey_spec != instrument_header:
			# 	self.survey_spec = instrument_header
			# 	print("[Spector] updating survey_spec to reflect instrument in spec.fits header -- {}".format(instrument_header))
				
			return at.Column(spec, name=['spec']), at.Column(ws, name=['ws']), at.Column(ivar, name=['ivar'])
		else: 
			raise NameError("[Spector] survey_spec not recognized")


	def _get_norm_trans_func(self, band='i'):
		"""
		return normalized transmission function and its wavelength coordinate
		the normalization is such that int{ trans{l} * dlnl} = 1.

		Params
		------
		self
		band='i'

		Return
		------
		trans  (array)
		ws_trans (array)
		"""
		trans, ws_trans = filters.filtertools.getNormTransFunc(band=band, survey=self.survey)
		return trans, ws_trans*u.AA


	def _get_norm_trans(self, wavelength, band='i', bounds_error=False):
		"""
		return normalized transmission function at a specific wavelenth, see filters.getNormTrans()
		the normalization is such that int{ trans{l} * dlnl} = 1.

		Params
		------
		self
		wave (float): wavelength to evaluate the function at
		band='i'
		bounds_error (bool): 
			whether to raise error when interpolation outside of array is attempted

		Return
		------
		trans  (float)
		"""
		trans = filters.filtertools.getNormTrans(wavelength, band=band, survey=self.survey, bounds_error=bounds_error) 
		return trans


	def _calc_Fnu_in_band(self, band, component='all'):
		"""
		Params
		------
		band=band
		component='all': from ['all', 'cont', 'line']
			which spectral component to operate on

		Return
		------
		Fnu (quantity in units "erg s-1 cm-2 Hz-1")
		"""

		spec, ws = self.get_spec_ws_from_spectab(component=component)

		trans, ws_trans = self._get_norm_trans_func(band=band)

		Fnu = filters.inttools.calc_Fnu_in_band_from_fl(fl=spec, ws=ws, trans=trans, ws_trans=ws_trans, isnormed=True)
		return Fnu


	def _calc_mAB_in_band(self, band, component='all'):
		Fnu = self._calc_Fnu_in_band(band=band, component='all')
		return Fnu.to(u.ABmag)


	def __get_specmag_colname(self, band, component, fluxquantity):
		"""
		band : e.g. 'i'
		component: from ['all', 'cont', 'line']
		fluxquantity: from ['fnu', 'mag']
		"""
		tag_component = {'all': '', 'cont': 'cont', 'line': 'line', 'contextrp': 'contextrp'}
		tag_quantity = {'fnu': 'Fnu', 'mag': 'Mag'}
		return 'spec{0}{1}_{2}'.format(tag_component[component], tag_quantity[fluxquantity], band)


	def __get_spectab_colname(self, component):
		"""
		band : e.g. 'i'
		component: from ['all', 'cont', 'line']
		fluxquantity: from ['fnu', 'mag']
		"""
		if component == 'ws':
			return 'ws'
		else: 
			tag_component = {'all': '', 'cont': 'cont', 'line': 'line', 'contextrp': 'contextrp'}
			return 'spec{0}'.format(tag_component[component])
