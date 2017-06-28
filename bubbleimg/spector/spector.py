# spector.py 
# ALS 2017/06/01

import numpy as np
import astropy.table as at
import matplotlib.pyplot as plt
import astropy.units as u

import os

from ..obsobj import Operator
from .. import filters
import getconti
# import inttools
import extrap
import linelist
import linefrac

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
			self.obj.add_sdss(toload_photoobj=False)
			if self.survey_spec == 'auto':
				self.survey_spec = self.obj.sdss.instrument.lower()

			# sanity check - spec name consistent
			if (self.survey_spec != self.obj.sdss.instrument.lower()):
				raise Exception("[spector] survey_spec inconsistent with sdss_xid.csv:instrument")

		# set z
		if hasattr(self.obj, 'z'):
			self.z = kwargs.pop('z', self.obj.z)
		elif self.survey_spec in ['sdss', 'boss', 'boss', 'auto']:
			self.obj.add_sdss(toload_photoobj=False)
			self.z = kwargs.pop('z', self.obj.sdss.z) 
		else: 
			self.z = kwargs.pop('z') 

		# set others
		self.bands = filters.filtertools.surveybands[self.survey]
		self.waverange = filters.filtertools.waverange[self.survey]

		# define paths
		self.fp_spec_decomposed = self.dir_obj+'spec_decomposed.ecsv'
		self.fp_spec_contextrp = self.dir_obj+'spec_contextrp.ecsv'
		self.fp_spec_mag = self.dir_obj+'spec_mag.csv'


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
		u_lcoord = u.AA
		"""

		fn = self.fp_spec_decomposed

		if not os.path.isfile(fn) or forceload_from_fits:
			return self.__read_spec_ws_from_fits()
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
		u_lcoord = u.AA
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
		"""

		fn = self.fp_spec_decomposed

		if (not os.path.isfile(fn)) or overwrite:
			spec, ws = self.get_spec_ws()
			speccon, specline, ws = getconti.decompose_cont_line_t2AGN(spec, ws, self.z)

			tab = at.Table([ws, spec, speccon, specline], names=['ws', 'spec', 'speccont', 'specline'])
			tab.write(fn, format='ascii.ecsv', overwrite=overwrite)

			# sanity check: units are identical
			units = [tab[col].unit for col in ['spec', 'speccont', 'specline']]
			if len(set(units)) > 1:
				raise Exception("[spector] units in table spec_decomposed are not identical")

		status = os.path.isfile(fn)
		return status 


	def make_spec_contextrp_ecsv(self, overwrite=False):
		""" extrapolate continuum to cover all of the wavelength range of filters """
		fn = self.fp_spec_contextrp

		if (not os.path.isfile(fn)) or overwrite:
			speccon, ws = self.get_spec_ws_from_spectab(component='cont')

			l0, l1 = self.waverange
			
			speccon_ext, lcoord_ext = extrap.extrap_to_end(ys=speccon, xs=ws, x_end=l0, polydeg=1, extbase_length=2000.)
			speccon_ext, lcoord_ext = extrap.extrap_to_end(ys=speccon_ext, xs=lcoord_ext, x_end=l1, polydeg=1, extbase_length=2000.)

			col = self.__get_spectab_colname('contextrp')
			tab = at.Table([lcoord_ext, speccon_ext], names=['ws', col])
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
			tabmag = at.Table()
			tabfnu = at.Table()

			for component in ['all', 'cont', 'line', 'contextrp']:
				for band in self.bands:
					colfnu = self.__get_specmag_colname(band, component=component, fluxquantity='fnu')
					colmag = self.__get_specmag_colname(band, component=component, fluxquantity='mag')
					try:
						fnu = self._calc_Fnu_in_band(band=band, component=component)
					except:
						print("[spector] skip calculating fnu of {} in band {}".format(component, band))
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
		calculate the ratio flux_line / fnu_band, for the conversion between band image and line image. 

		by multiplying this ratio to the band image (in fnu [erg/s/cm^2/Hz]) one can obtain the observed flux of the line (in f [erg/s/cm^2]). Notice that this is flux in observed frame, for measurements one still needs to transform it to the rest frame (depending on z). 

		Params
		------
		self
		band (str)
		line (str)

		Return
		------
		ratio (quantity in unit of Hz)
		"""
		print("WARNING: to write error propogation as well.......")

		fl, flerr = self._get_line_flux(line=line, wunit=False)
		wl = self._get_line_obs_wave(line=line, wunit=False)
		Tl = self._get_norm_trans(wavelength=wl, band=band, bounds_error=True)

		# assemble line list
		lines = self._list_stronglines_in_band(band=band)
		fs = np.ndarray(len(lines))
		ws = np.ndarray(len(lines))
		Ts = np.ndarray(len(lines))
		for i, l in enumerate(lines):
			f, ferr = self._get_line_flux(line=l, wunit=False)
			w = self._get_line_obs_wave(line=l, wunit=False)
			T = self._get_norm_trans(wavelength=w, band=band, bounds_error=True)
			fs[i] = f
			ws[i] = w
			Ts[i] = T

		ratio = linefrac.fline_over_fnuband(fl, wl, Tl, fs, ws, Ts)

		ratio = ratio.to(u.Hz)

		return ratio


	def _get_line_flux(self, line='OIII5008', wunit=False): 
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


	def __read_spec_ws_from_fits(self):
		"""
		Return
		------
		spec (nparray)
		ws (nparray)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_lcoord = u.AA
		"""
		if self.survey_spec in ['sdss', 'boss', 'boss']:
			spec, ws = self.obj.sdss.get_speclcoord(wunit=True)
			return at.Column(spec, name=['spec']), at.Column(ws, name=['ws'])
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
