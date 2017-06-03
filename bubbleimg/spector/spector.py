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

		z (float):  
			redshift, if not provided, use obj.sdss.z
		survey = 'hsc' (str): 
			survey of the photometric system
			if not provided, use self.obj.survey (or 'hsc' if attribute does not exist)
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
		"""
		
		super(Spector, self).__init__(**kwargs)

		if hasattr(self.obj, 'survey'):
			default_survey = self.obj.survey
		else: 
			default_survey = 'hsc'

		self.survey = kwargs.pop('survey', default_survey)
		self.survey_spec = kwargs.pop('survey_spec', 'auto')
		self.bands = filters.filtertools.surveybands[self.survey]

		if self.survey_spec in ['sdss', 'boss', 'boss', 'auto']:
			self.obj.add_sdss(toload_photoobj=False)
			self.z = kwargs.pop('z', self.obj.sdss.z) 

			if self.survey_spec == 'auto':
				self.survey_spec = self.obj.sdss.instrument.lower()

			# sanity check - spec name consistent
			if (self.survey_spec != self.obj.sdss.instrument.lower()):
				raise Exception("[spector] survey_spec in consistent with sdss_xid.csv:instrument")

		else: 
			self.z = kwargs.pop('z') 

		# define paths
		self.fp_spec_decomposed = self.dir_obj+'spec_decomposed.ecsv'
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
			speccon, specline, lcoord = getconti.decompose_conti_line_t2AGN(spec, lcoord, self.z)

			tab = at.Table([lcoord, spec, speccon, specline], names=['lcoord', 'spec', 'specconti', 'specline'])
			tab.write(fn, format='ascii.ecsv', overwrite=overwrite)

			# sanity check: units are identical
			units = [tab[col].unit for col in ['spec', 'specconti', 'specline']]
			if len(set(units)) > 1:
				raise Exception("[spector] units in table spec_decomposed are not identical")

		status = os.path.isfile(fn)
		return status 


	def get_specconti_lcoord(self):
		"""
		Return
		------
		specconti (astropy table column with unit)
		lcoord (astropy table column with unit)

		Default units
		-------------
		u_spec = 1.e-17*u.Unit('erg / (Angstrom cm2 s)')
		u_lcoord = u.AA
		"""
		fn = self.fp_spec_decomposed
		self.make_spec_decomposed_ecsv(overwrite=False)

		tab = at.Table.read(fn, format='ascii.ecsv')

		return tab['specconti'], tab['lcoord']


	def get_specline_lcoord(self):
		"""
		Return
		------
		specline (astropy table column with unit)
		lcoord (astropy table column with unit)
		"""
		fn = self.fp_spec_decomposed
		self.make_spec_decomposed_ecsv(overwrite=False)

		tab = at.Table.read(fn, format='ascii.ecsv')

		return tab['specline'], tab['lcoord']
	

	def get_Fnu_in_band(self, band, spec_component='all'):
		"""
		Params
		------
		band=band
		spec_component='all': from ['all', 'conti', 'line']
			which spectral component to operate on

		Return
		------
		Fnu (quantity in units "erg s-1 cm-2 Hz-1")
		"""

		if spec_component == 'all':
			spec, lcoord = self.get_spec_lcoord()
		elif spec_component == 'conti':
			spec, lcoord = self.get_specconti_lcoord()
		elif spec_component == 'line':
			spec, lcoord = self.get_specline_lcoord()
		else:
			raise Exception("[spector] spec_component not recognized")

		trans, lcoord_trans = self.get_filter_trans(band=band)

		Fnu = inttools.calc_Fnu_in_band_from_fl(fl=spec, ls=lcoord, trans=trans, ls_trans=lcoord_trans)

		return Fnu


	def get_mAB_in_band(self, band, spec_component='all'):
		Fnu = self.get_Fnu_in_band(band=band, spec_component='all')
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
		def get_colname(band, spec_component, stype):
			"""
			band : e.g. 'i'
			spec_component: from ['all', 'conti', 'line']
			stype: from ['Fnu', 'Mag']
			"""
			scomponent = {'all': '', 'conti': 'conti', 'line': 'line'}
			return 'spec{0}{1}_{2}'.format(scomponent[spec_component], stype, band)

		#==========================================================================

		fn = self.fp_spec_mag

		if not os.path.isfile(fn) or overwrite:
			tabmag = at.Table()
			tabfnu = at.Table()

			for spec_component in ['all', 'conti', 'line']:
				for band in self.bands:
					colfnu = get_colname(band, spec_component, stype='Fnu')
					colmag = get_colname(band, spec_component, stype='Mag')
					fnu = self.get_Fnu_in_band(band=band, spec_component=spec_component)
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


	def plot_spec(self, wfilters=True, wconti=False, overwrite=False):
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
			plt.figure(1, figsize=(12, 6))
			plt.clf()
			if wfilters:
				for band in self.bands:
					trans, lcoord_trans = self.get_filter_trans(band=band)
					plt.plot(lcoord_trans, trans/max(trans), label=band)

			spec, lcoord = self.get_spec_lcoord()
			norm = max(spec)
			
			plt.plot(lcoord, spec/norm, color='black', lw=2, label='__nolabel__')
			plt.plot(0., 0., ls='', color='black', label='z='+'%.3f'%self.z) # show z in legend

			if wconti:
				specconti, lcoord = self.get_specconti_lcoord()
				plt.plot(lcoord, specconti/norm, color='blue', lw=2, label='__nolabel__')

			plt.legend(loc='upper right')
			if self.survey_spec == 'sdss':
				plt.xlim(2980.0, 11230.0)
			elif self.survey_spec == 'boss':
				plt.xlim(3500.0, 12000.0)
			plt.ylim(0., 1.)
			plt.savefig(fn)

		status = os.path.isfile(fn)
		return status 



