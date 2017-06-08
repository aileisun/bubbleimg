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
import inttools
import extrap

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
				raise Exception("[spector] survey_spec in consistent with sdss_xid.csv:instrument")

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


	def __read_spec_lcoord_from_fits(self):
		"""
		Return
		------
		spec (nparray)
		lcoord (nparray)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_lcoord = u.AA
		"""
		if self.survey_spec in ['sdss', 'boss', 'boss']:
			spec, lcoord = self.obj.sdss.get_speclcoord(wunit=True)
			return at.Column(spec, name=['spec']), at.Column(lcoord, name=['lcoord'])
		else: 
			raise NameError("[Spector] survey_spec not recognized")


	def get_spec_lcoord(self, forceload_from_fits=False):
		"""
		read spec and lcoord, either from spec_decomposed.ecsv or spec.fits (if forced or spec_decomposed.ecsv does not exist)

		Param
		-----
		forceload_from_fits=False:
			if true, then load from fits. 

		Return
		------
		spec (astropy table column with unit)
		lcoord (astropy table column with unit)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_lcoord = u.AA
		"""

		fn = self.fp_spec_decomposed

		if not os.path.isfile(fn) or forceload_from_fits:
			return self.__read_spec_lcoord_from_fits()
		else: 
			tab = at.Table.read(fn, format='ascii.ecsv')
			return tab['spec'], tab['lcoord']


	def get_spec_lcoord_from_spectab(self, component, fn=None):
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
		lcoord (astropy table column with unit)

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

		return tab[col], tab['lcoord']


	def get_filter_trans(self, band='i'):
		"""
		Params
		------
		self
		band='i'

		Return
		------
		trans  (array)
		lcoord_trans (array)
		"""
		trans, lcoord_trans = filters.getFilterResponseFunc(band=band, survey=self.survey)
		return trans, lcoord_trans*u.AA


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
			spec, lcoord = self.get_spec_lcoord()
			speccon, specline, lcoord = getconti.decompose_cont_line_t2AGN(spec, lcoord, self.z)

			tab = at.Table([lcoord, spec, speccon, specline], names=['lcoord', 'spec', 'speccont', 'specline'])
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
			speccon, lcoord = self.get_spec_lcoord_from_spectab(component='cont')

			l0, l1 = self.waverange
			
			speccon_ext, lcoord_ext = extrap.extrap_to_end(ys=speccon, xs=lcoord, x_end=l0, polydeg=1, extbase_length=1000.)
			speccon_ext, lcoord_ext = extrap.extrap_to_end(ys=speccon_ext, xs=lcoord_ext, x_end=l1, polydeg=1, extbase_length=1000.)

			col = self.__get_spectab_colname('contextrp')
			tab = at.Table([lcoord_ext, speccon_ext], names=['lcoord', col])
			tab.write(fn, format='ascii.ecsv', overwrite=overwrite)

		status = os.path.isfile(fn)
		return status 


	def calc_Fnu_in_band(self, band, component='all'):
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

		spec, lcoord = self.get_spec_lcoord_from_spectab(component=component)

		trans, lcoord_trans = self.get_filter_trans(band=band)

		Fnu = inttools.calc_Fnu_in_band_from_fl(fl=spec, ls=lcoord, trans=trans, ls_trans=lcoord_trans)

		return Fnu


	def calc_mAB_in_band(self, band, component='all'):
		Fnu = self.calc_Fnu_in_band(band=band, component='all')
		return Fnu.to(u.ABmag)


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
						fnu = self.calc_Fnu_in_band(band=band, component=component)
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
		if component == 'lcoord':
			return 'lcoord'
		else: 
			tag_component = {'all': '', 'cont': 'cont', 'line': 'line', 'contextrp': 'contextrp'}
			return 'spec{0}'.format(tag_component[component])


	def get_spec_mag_tab(self):
		""" return spec_mag table"""
		self.make_spec_mag(overwrite=False)
		return at.Table.read(self.fp_spec_mag, comment='#', format='ascii.csv')


	def get_spec_mag_value(self, band, component='all', fluxquantity='mag'):
		""" 
		return requested value, reading from spec_mag.csv

		Params
		------
		component='line'
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



	def plot_spec(self, wfilters=True, wconti=False, wcontextrp=False, wline=False, overwrite=True):
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
			if wfilters:
				for band in self.bands:
					trans, lcoord_trans = self.get_filter_trans(band=band)
					plt.plot(lcoord_trans, trans/max(trans), label=band)

			spec, lcoord = self.get_spec_lcoord()
			norm = max(spec)
			
			plt.plot(lcoord, spec/norm, color='0.3', lw=1.5, label='__nolabel__')
			plt.plot(0., 0., ls='', color='black', label='z='+'%.3f'%self.z) # show z in legend

			if wline:
				speccont, lcoord = self.get_spec_lcoord_from_spectab(component='cont')
				specline, lcoord = self.get_spec_lcoord_from_spectab(component='line')
				specplot = specline+speccont
				specplot[specline==0] = np.nan
				plt.plot(lcoord, specplot/norm, color='black', lw=1.5, label='line')

			if wconti:
				speccont, lcoord = self.get_spec_lcoord_from_spectab(component='cont')
				plt.plot(lcoord, speccont/norm, color='cyan', lw=2, label='continuum')

			if wcontextrp:
				speccextrp, lcoord = self.get_spec_lcoord_from_spectab(component='contextrp')
				plt.plot(lcoord, speccextrp/norm, color='grey', lw=0.5, label='extrapolated conti')


			plt.legend(loc='upper right')
			if self.survey_spec == 'sdss':
				plt.xlim(2980.0, 11230.0)
			elif self.survey_spec == 'boss':
				plt.xlim(3500.0, 12000.0)
			plt.ylim(0., 1.)
			plt.savefig(fn, overwrite=overwrite)

		status = os.path.isfile(fn)
		return status 




	# def get_speccont_lcoord(self):
	# 	"""
	# 	Return
	# 	------
	# 	speccont (astropy table column with unit)
	# 	lcoord (astropy table column with unit)

	# 	Default units
	# 	-------------
	# 	u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
	# 	u_lcoord = u.AA
	# 	"""
	# 	fn = self.fp_spec_decomposed
	# 	self.make_spec_decomposed_ecsv(overwrite=False)

	# 	tab = at.Table.read(fn, format='ascii.ecsv')

	# 	return tab['speccont'], tab['lcoord']


	# def get_specline_lcoord(self):
	# 	"""
	# 	Return
	# 	------
	# 	specline (astropy table column with unit)
	# 	lcoord (astropy table column with unit)
	# 	"""
	# 	fn = self.fp_spec_decomposed
	# 	self.make_spec_decomposed_ecsv(overwrite=False)

	# 	tab = at.Table.read(fn, format='ascii.ecsv')

	# 	return tab['specline'], tab['lcoord']
	
